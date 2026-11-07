.PHONY: init plan test

CONFIG ?= config/puget_sound.toml

init:
	python3 -m puget_modflow init-config --output $(CONFIG)

plan:
	python3 -m puget_modflow plan --config $(CONFIG) --output runs/puget_sound/pipeline_plan.json

test:
	PYTHONPATH=src python3 -m unittest discover -s tests
