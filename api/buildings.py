# Vercel 서버 함수: 라우팅·계산은 루트 viewimpact.py 의 Handler 가 담당
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from viewimpact import Handler as handler  # noqa: E402,F401
