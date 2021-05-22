# ViennaRNA container must be built first, so we remove it from the list; then prepend it.
DOCKER_BUILDS := vienna_rna $(filter-out vienna_rna,$(patsubst dockers/%/,%, $(dir $(wildcard dockers/*/Dockerfile))))

all: create_local_config fetch_python_requirements build_dockers

create_local_config:
	@echo "CREATING LOCAL CONFIG FILE"
	cp config.json -n config_local.json  # Shouldn't overwrite existing local config

fetch_python_requirements:
	@echo "FETCHING PYTHON REQUIREMENTS"
	pip install -r requirements.txt

build_docker_%:
	@echo "BUILDING IMAGE '$*'"
	docker build -t $* --force-rm dockers/$*/

build_dockers: $(addprefix build_docker_,$(DOCKER_BUILDS))
