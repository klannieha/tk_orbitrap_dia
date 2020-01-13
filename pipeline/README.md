Pipeline of running DIA analysis:

1. Convert raw files to mzML using ProteoWizard MSConvert (locally or docker)
	- script: ./conversion.sh
	- run script: ./run_conversion.sh

2. Create spectral library from MQ DDA fractionated data
	- script: ./create_library.py

3. Generate OSW spectral library with commands:
	OpenSwathAssayGenerator --in mqout.tsv --out spectral_lib.pqp 
	OpenSwathDecoyGenerator --in spectral_lib.pqp --out spectral_lib_target_decoy.pqp --switchKR true --method pseudo-reverse

4. Run OpenSwathWorkflow:
	- script: ./run_osw_HEK.sh
	- run script: ./run_osw.sh

5. Run statistical analysis with PyProphet 
	- script: ./run_pyprophet.sh

