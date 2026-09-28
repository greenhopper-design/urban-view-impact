"""Urban View Impact — 신축 전·후 주변 건물의 층별 조망률 변화 분석.

    python viewimpact.py demo    # 합성 데이터로 자체 검증 (API 키 불필요)
    python viewimpact.py serve   # http://localhost:8000 지도 앱

환경변수
    VWORLD_KEY     브이월드 인증키 (serve 에 필요)
    VWORLD_DOMAIN  인증키 발급 시 등록한 도메인 (기본 localhost)
    VIEW_DEM       DEM GeoTIFF 경로 (선택, 없으면 모든 지반고 0m = 평지 가정)

높이 모델 (모두 해발고도)
    관찰점  = 건물 GL + (층-1) x 층고 + 1.5m
    건물지붕 = 건물 GL + 지상층수 x 층고
    조망대상 = 대상 GL + 대상 높이
"""
import json
import math
import os
import sys
from functools import lru_cache
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import numpy as np
import requests
import shapely
from shapely.geometry import LineString, Point, shape
from shapely.ops import nearest_points
from shapely.strtree import STRtree

EYE = 1.5            # 바닥에서 눈높이 (m)
EARTH_R = 6_371_000  # 지구 곡률 보정용
TERRAIN_STEP = 30    # 지형 가림 검사 간격 (m)
TILE = 0.01          # 브이월드 조회 타일 (도), 약 0.9 x 1.1 km < 2 km² 제한
VWORLD = "https://api.vworld.kr/req/data"
SEOUL = (126.76, 37.41, 127.19, 37.72)  # 분석 범위 (경도·위도 min/max)
HERE = Path(__file__).parent


# ---------------------------------------------------------------- 좌표
class Local:
    """위경도 -> 부지 중심 기준 로컬 미터 좌표.
    ponytail: 등장방형 근사, 수 km 이내 오차 0.1% 미만. 더 넓으면 pyproj EPSG:5186 사용."""

    def __init__(self, lon0, lat0):
        self.o = np.array([lon0, lat0])
        self.k = np.array([111_320 * math.cos(math.radians(lat0)), 110_574])

    def geom(self, g):
        return shapely.transform(g, lambda c: (c - self.o) * self.k)

    def xy(self, lon, lat):
        return tuple((np.array([lon, lat]) - self.o) * self.k)

    def lonlat(self, xy):
        return np.asarray(xy) / self.k + self.o


# ---------------------------------------------------------------- 데이터
def fetch_buildings(bbox):
    """bbox(경도,위도 min/max) 안의 브이월드 도로명주소건물(LT_C_SPBD). 서울 밖은 잘라냄."""
    minx, miny = max(bbox[0], SEOUL[0]), max(bbox[1], SEOUL[1])
    maxx, maxy = min(bbox[2], SEOUL[2]), min(bbox[3], SEOUL[3])
    if minx >= maxx or miny >= maxy:
        return []
    ix = range(math.floor(minx / TILE), math.floor(maxx / TILE) + 1)
    iy = range(math.floor(miny / TILE), math.floor(maxy / TILE) + 1)
    if len(ix) * len(iy) > 30:
        raise ValueError("조회 범위가 너무 넓습니다 (조망대상을 더 가까이 두세요)")
    out = {}
    for i in ix:
        for j in iy:
            for f in _fetch_tile(i, j):  # 타일 경계 건물 중복 제거
                out[f["id"]] = f
    return list(out.values())


@lru_cache(maxsize=512)
def _fetch_tile(i, j):
    box = (i * TILE, j * TILE, (i + 1) * TILE, (j + 1) * TILE)
    out, page = [], 1
    while True:
        r = requests.get(VWORLD, timeout=30, params=dict(
            service="data", request="GetFeature", data="LT_C_SPBD",
            key=os.environ["VWORLD_KEY"], domain=os.environ.get("VWORLD_DOMAIN", "localhost"),
            geomFilter="BOX(%.6f,%.6f,%.6f,%.6f)" % box, crs="EPSG:4326",
            size=1000, page=page, format="json", geometry="true", attribute="true",
        )).json()["response"]
        if r["status"] == "NOT_FOUND":
            return tuple(out)
        if r["status"] != "OK":
            raise RuntimeError("브이월드 오류: %s" % r.get("error", r["status"]))
        for f in r["result"]["featureCollection"]["features"]:
            p = f["properties"]
            out.append({
                "id": str(f.get("id") or p.get("bd_mgt_sn")),
                "name": " ".join(filter(None, [p.get("buld_nm"), p.get("buld_nm_dc")])),
                "floors": int(float(p.get("gro_flo_co") or 0)),
                "geometry": f["geometry"],
            })
        if page >= int(r["page"]["total"]):
            return tuple(out)
        page += 1


class Ground:
    """DEM에서 지반고 조회. DEM이 없으면 0 (평지 가정)."""

    def __init__(self, path=None):
        self.ds = None
        if path:
            import rasterio
            from rasterio.warp import transform
            self.ds, self._tr = rasterio.open(path), transform

    def at(self, lonlats):
        lonlats = np.asarray(lonlats, float).reshape(-1, 2)
        if self.ds is None:
            return np.zeros(len(lonlats))
        xs, ys = self._tr("EPSG:4326", self.ds.crs, lonlats[:, 0].tolist(), lonlats[:, 1].tolist())
        v = np.array([s[0] for s in self.ds.sample(zip(xs, ys))], float)
        if self.ds.nodata is not None:
            v[v == self.ds.nodata] = np.nan
        return np.nan_to_num(v, nan=0.0)

    def gl(self, geom_lonlat):
        """건물 GL = 외곽선 꼭짓점 + 중심의 DEM 평균 (경사지 평균 지표면 근사)."""
        pts = np.vstack([shapely.get_coordinates(geom_lonlat), [geom_lonlat.centroid.coords[0]]])
        return float(self.at(pts).mean())


# ---------------------------------------------------------------- 분석
def _span(seg, geom, p):
    """선분이 도형을 지나는 구간의 P로부터 거리 (최소, 최대). 안 지나면 None."""
    c = shapely.get_coordinates(seg.intersection(geom))
    if not len(c):
        return None
    d = np.hypot(c[:, 0] - p[0], c[:, 1] - p[1])
    return d.min(), d.max()


def _blocked(zp, zt, D, hits):
    """hits: [(거리, 가림높이)] 중 하나라도 시선보다 높으면 가려짐. 지구 곡률 반영."""
    for d, z in hits:
        line = zp + (zt - zp) * d / D
        if z - d * (D - d) / (2 * EARTH_R) > line:
            return True
    return False


def analyze(buildings, site_id, new_height, targets, floor_h=3.0, radius=300, ground=None, gl=None):
    """
    buildings: [{id, name, floors, geometry(GeoJSON, 경위도)}]
    targets:   [{lon, lat, h, gl?}]  h = 대상 지반 위 높이 (한강 0, 남산타워 236 등)
    gl:        {건물id: 지반고} 브라우저가 브이월드 지형에서 뽑은 값. 없으면 DEM/0
    반환: 반경 내 주변 건물별 층별 조망률(%) 신축 전/후, 최대 감소폭 순 정렬
    """
    if not targets:
        raise ValueError("조망대상을 하나 이상 지정하세요")
    ground = ground or Ground()
    ids = [b["id"] for b in buildings]
    if site_id not in ids:
        raise ValueError("신축 부지 건물을 찾을 수 없습니다")
    s = ids.index(site_id)

    ll = [shape(b["geometry"]) for b in buildings]
    L = Local(*ll[s].centroid.coords[0])
    polys = [L.geom(g) for g in ll]
    gl = [gl[b["id"]] if gl and b["id"] in gl else ground.gl(g) for b, g in zip(buildings, ll)]
    floors = [max(b["floors"], 1) for b in buildings]  # 층수 누락 건물은 1층 가정
    roof = [gl[i] + floors[i] * floor_h for i in range(len(buildings))]
    old_site, new_site = roof[s], gl[s] + new_height

    T = [(np.array(L.xy(t["lon"], t["lat"])),
          (t["gl"] if "gl" in t else float(ground.at([[t["lon"], t["lat"]]])[0])) + t.get("h", 0))
         for t in targets]
    tree = STRtree(polys)
    site = polys[s]

    results = []
    for i, poly in enumerate(polys):
        if i == s or poly.distance(site) > radius:
            continue
        # (관찰건물, 대상)별 가림 요소는 층과 무관 -> 한 번만 계산
        rays = []
        for txy, zt in T:
            p = np.array(nearest_points(poly.boundary, Point(txy))[0].coords[0])
            u = txy - p
            D = float(np.hypot(*u))
            if D < 1:
                continue
            p = p + u / D * 0.5  # 자기 벽 바깥으로 0.5m
            seg = LineString([p, txy])
            others = []
            for j in tree.query(seg, predicate="intersects"):
                if j in (i, s):
                    continue
                others += [(d, roof[j]) for d in _span(seg, polys[j], p) or ()]
            if ground.ds is not None:  # 중간 지형(언덕) 가림. ponytail: 브이월드 지형 GL만 쓰면 언덕 가림은 생략
                ds = np.arange(TERRAIN_STEP, D - TERRAIN_STEP, TERRAIN_STEP)[:200]
                if len(ds):
                    zs = ground.at(L.lonlat(p + np.outer(ds / D, u)))
                    others += list(zip(ds, zs))
            sp = _span(seg, site, p)
            rays.append((zt, D, others, sp or ()))

        rows = []
        for f in range(1, floors[i] + 1):
            zp = gl[i] + (f - 1) * floor_h + EYE
            before = after = 0
            for zt, D, others, sp in rays:
                if _blocked(zp, zt, D, others):
                    continue
                before += not _blocked(zp, zt, D, [(d, old_site) for d in sp])
                after += not _blocked(zp, zt, D, [(d, new_site) for d in sp])
            n = len(T)
            rows.append({"floor": f, "before": round(100 * before / n, 1),
                         "after": round(100 * after / n, 1),
                         "change": round(100 * (after - before) / n, 1)})
        results.append({"id": ids[i], "name": buildings[i]["name"], "gl": round(gl[i], 1),
                        "max_drop": -min(r["change"] for r in rows), "floors": rows})
    results.sort(key=lambda r: -r["max_drop"])
    return results


# ---------------------------------------------------------------- 서버
class Handler(BaseHTTPRequestHandler):
    ground = Ground()  # serve() 에서 DEM 지정 시 교체

    def _send(self, code, body, ctype="application/json; charset=utf-8"):
        body = body if isinstance(body, bytes) else json.dumps(body, ensure_ascii=False).encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        u = urlparse(self.path)
        try:
            if u.path == "/":
                # 브이월드 3D 스크립트는 브라우저에서 키가 필요 (도메인 제한 키라 노출 허용 범위)
                html = (HERE / "index.html").read_text(encoding="utf-8")
                html = html.replace("{{VWORLD_KEY}}", os.environ["VWORLD_KEY"]).replace(
                    "{{VWORLD_DOMAIN}}", os.environ.get("VWORLD_DOMAIN", "localhost"))
                return self._send(200, html.encode(), "text/html; charset=utf-8")
            if u.path == "/api/buildings":
                bbox = tuple(map(float, parse_qs(u.query)["bbox"][0].split(",")))
                return self._send(200, {"buildings": fetch_buildings(bbox), "dem": self.ground.ds is not None})
            self._send(404, {"error": "not found"})
        except Exception as e:
            self._send(500, {"error": str(e)})

    def do_POST(self):
        try:
            q = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            res = analyze(fetch_buildings(q["bbox"]), q["site_id"], float(q["new_height"]), q["targets"],
                          float(q["floor_h"]), float(q["radius"]), self.ground, q.get("gl"))
            self._send(200, {"results": res})
        except Exception as e:
            self._send(400, {"error": str(e)})


def serve(port=8000):
    if not os.environ.get("VWORLD_KEY"):
        sys.exit("VWORLD_KEY 환경변수를 설정하세요")
    Handler.ground = Ground(os.environ.get("VIEW_DEM"))
    print("http://localhost:%d  (DEM: %s)" % (port, os.environ.get("VIEW_DEM") or "없음, 평지 가정"))
    ThreadingHTTPServer(("127.0.0.1", port), Handler).serve_forever()


# ---------------------------------------------------------------- 자체 검증
def demo():
    """관찰건물 A(10층) -- 부지 S -- 조망대상 T(한강, 남쪽 300m) 일렬 배치."""
    lon0, lat0 = 127.0, 37.53
    kx, ky = 111_320 * math.cos(math.radians(lat0)), 110_574

    def box(x0, y0, x1, y1):
        return {"type": "Polygon", "coordinates": [[(lon0 + x / kx, lat0 + y / ky) for x, y in
                                                    [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]]]}

    bs = [{"id": "A", "name": "관찰동", "floors": 10, "geometry": box(0, 0, 20, 10)},
          {"id": "S", "name": "부지", "floors": 2, "geometry": box(0, 50, 20, 60)}]
    t = [{"lon": lon0 + 10 / kx, "lat": lat0 + 300 / ky, "h": 0}]

    r = analyze(bs, "S", 60, t)[0]["floors"]  # 2층(6m) -> 60m
    assert r[0]["before"] == 0, r[0]           # 1층: 기존 2층 건물에도 이미 가림
    assert r[9]["before"] == 100, r[9]         # 10층: 기존엔 보임
    assert r[9]["after"] == 0, r[9]            # 10층: 신축 60m 후 가림
    assert analyze(bs, "S", 6, t)[0]["max_drop"] == 0  # 높이 그대로면 변화 없음
    hill = analyze(bs, "S", 60, t, gl={"A": 60, "S": 0})[0]["floors"]
    assert hill[9]["after"] == 100, hill[9]    # 관찰동이 60m 언덕 위면 신축 60m 넘어 보임
    for row in r:
        print("%2dF  기존 %5.1f%%  신축후 %5.1f%%  변화 %+.1f%%p" % (row["floor"], row["before"], row["after"], row["change"]))
    print("demo OK")


if __name__ == "__main__":
    {"demo": demo, "serve": serve}.get(sys.argv[1] if len(sys.argv) > 1 else "", lambda: sys.exit(__doc__))()
