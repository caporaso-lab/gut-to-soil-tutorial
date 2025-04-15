# helpers
_copy-env-file:
	cp environment-files/readthedocs.yml docs/_static/environment.yml

_copy-data:
	@if [ -d "docs/data/" ]; then \
		cp -r docs/data/ docs/_build/html/data/; \
		echo "Copied data directory."; \
	else \
		echo "docs/data/ not found, skipping copy."; \
	fi

_build-html:
	cd docs && jupyter book build --html

_build-fast-preview:
	cd docs && Q2DOC_FASTMODE= jupyter book build --html

_build-preview:
	cd docs && Q2DOC_PREVIEW= jupyter book build --html

# main targets
html: _copy-env-file _build-html _copy-data

fast-preview: _copy-env-file _build-fast-preview _copy-data

preview: _copy-env-file _build-preview _copy-data

serve:
	npx serve docs/_build/html/ -p 4000

clean:
	rm -rf docs/_build/html/