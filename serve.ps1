# Discern Quick Launcher (Windows PowerShell)
Write-Host "Starting Discern at http://127.0.0.1:8000 ..." -ForegroundColor Cyan
python -m uvicorn discern.api.app:app --app-dir src --host 127.0.0.1 --port 8000 --no-access-log
