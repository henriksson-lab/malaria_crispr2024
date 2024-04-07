docker:
	docker build -t malaria_crispr2024 .
docker_testrun:
	docker run --rm -p 3838:3838 malaria_crispr2024
