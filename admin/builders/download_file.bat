:: download a file using curl or PowerShell
SET URL=%~1
SET OUT=%~2

echo Downloading %URL% to %OUT%
where curl >nul 2>&1
IF %ERRORLEVEL%==0 (
    curl -L "%URL%" -o "%OUT%"
) ELSE (
    PowerShell "[Net.ServicePointManager]::SecurityProtocol = 'tls12, tls11, tls'; (new-object System.Net.WebClient).DownloadFile('%URL%','%OUT%')"
)