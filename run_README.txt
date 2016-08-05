
1) Setting up a new model

- Open code/wrapper.R and fill out model description
	* Specify user_options_file, stored in screentreat/code
	* Specify input_data_file, stored in screentreat/examples

- open user options file and set model options
	* change model_version in the top section to reflect 
		what you put in the wrapper
	* if making a small modification to another version, vimdiff this file
		with the user_options.R file from that one 
		(stored in examples/[version]/input) to make sure everything
		else is indeed the same

- Open input data file and fill out treatment info/make sure it's the right one

- Source wrapper.R to set up (run=FALSE) or set up and run (run=TRUE) the model



2) Re-running the model with changes

- Manually edit version_guide.csv (within screentreat) and model_info.csv (within the screentreat/examples/[model] folder)
- Rerun model using 1 of 2 options:

* OPTION 1 (preferred, I think)
- Edit the user_options.R file that is stored within screentreat/examples/[model]/input
- If necessary, edit the input.csv file stored in that same folder
- Then source that file to re-run model

* OPTION 2 
- Edit code/user_options.R, making sure the model_version variable at the top is correct
- Source code/wrapper.R to re-run the model

* DO NOT:
- Edit code/user_options.R and then source it. This doesn't update the code stored within the model folder, so they get out of sync.