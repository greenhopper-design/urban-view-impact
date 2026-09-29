# Vercel 서버 함수: 라우팅·계산은 루트 viewimpact.py 의 Handler 가 담당
# (Vercel 은 파일 안에 직접 정의된 handler 클래스만 함수로 인식)
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from viewimpact import Handler  # noqa: E402


class handler(Handler):
    pass
