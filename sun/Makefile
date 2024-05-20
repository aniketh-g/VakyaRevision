.PHONY: generate_data
generate_data:
	@python3 ./generate_data.py

.PHONY: revision
revision:
	@python3 ./siddhanta_revision.py

.PHONY: vakyas
vakyas:
	@python3 ./vakya_generator.py

.PHONY: historical_comparison
historical_comparison:
	@python3 ./siddhanta_historical_comparison.py

DATADIR = ./sun_computed_data

.PHONY: clean
clean:
	@rm -rf $(DATADIR)
	@mkdir $(DATADIR)
	@cd __pycache__; rm -rf *.pyc