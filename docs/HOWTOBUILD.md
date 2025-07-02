

- cd to `pressio-ops/docs`

- create a virtual env and activate it
```
python3 -m venv venv
```

- once activated, do:
```
pip install -r build_requirements.txt
```

- then build docs with:
```
sphinx-autobuild source/ _build/html
```
this will generated the html inside the `_build/html` directory and start a server.
Just open it and start watching.
You can now edit the source files inside the `docs/source` and sphinx will automatically update the website.

More info at: https://github.com/sphinx-doc/sphinx-autobuild
