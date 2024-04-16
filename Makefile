docker:
	docker build -t malaria_crispr2024 .
docker_testrun:
	docker run --rm -p 3838:3838 malaria_crispr2024
docker_push:
	#docker tag 045a4cba1107 mahogny83/malaria_crispr2024:20240416-095800
	#docker image push mahogny83/malaria_crispr2024:20240416-095800
	#need to update tag
