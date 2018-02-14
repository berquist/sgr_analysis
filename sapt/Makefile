extract:
	# This extracts the SAPT data.
	# python collect_sapt_data.py
	# No, *this* extracts the SAPT data.
	bash run.sh
	# This extracts the ALMO data. It also reads SAPT CSV files, so
	# run it after the SAPT parsing.
	python sapt_vs_almo.py
	python summary.py

plot:
	python cp_geometries.py
	python trellis.py
