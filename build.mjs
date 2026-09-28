// Vercel 빌드: index.html 의 {{VWORLD_KEY}} 를 환경변수로 치환해 public/ 에 출력
// (브이월드 3D 스크립트는 브라우저에서 키가 필요. 키는 등록 도메인에서만 동작)
import fs from 'node:fs'

const key = process.env.VWORLD_KEY
if (!key) throw new Error('VWORLD_KEY 환경변수가 없습니다')
fs.mkdirSync('public', { recursive: true })
fs.writeFileSync('public/index.html', fs.readFileSync('index.html', 'utf8').replaceAll('{{VWORLD_KEY}}', key))
