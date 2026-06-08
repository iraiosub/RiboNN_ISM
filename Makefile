install:
	mamba env create -f environment.yml -y

train:
	export MLFLOW_TRACKING_URI=sqlite:///mlruns.db
	python3 -m src.main --train

transfer_learning:
	@if [ ! -d "models" ]; then \
		echo "Downloading pretrained models..."; \
		mkdir -p tmp; \
		wget -t 0 -O tmp/weights.zip https://zenodo.org/records/17258709/files/weights.zip?download=1; \
		unzip -oq tmp/weights.zip -d models/; \
		rm -r tmp; \
	fi;
	export MLFLOW_TRACKING_URI=sqlite:///mlruns.db; \
	python3 -m src.main --transfer_learning

predict_human:
	@if [ ! -d "models/human" ]; then \
		echo "Downloading pretrained models..."; \
		mkdir -p tmp; \
		wget -t 0 -O tmp/weights.zip https://zenodo.org/records/17258709/files/weights.zip?download=1; \
		unzip -oq tmp/weights.zip -d models/; \
		rm -r tmp; \
	fi;
	python3 -m src.main --predict human

predict_mouse:
	@if [ ! -d "models/mouse" ]; then \
		echo "Downloading pretrained models..."; \
		mkdir -p tmp; \
		wget -t 0 -O tmp/weights.zip https://zenodo.org/records/17258709/files/weights.zip?download=1; \
		unzip -oq tmp/weights.zip -d models/; \
		rm -r tmp; \
	fi;
	python3 -m src.main --predict mouse
