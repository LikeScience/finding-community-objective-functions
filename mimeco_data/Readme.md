# Recovering objectives in _L. plantarum & A. muciniphila community_

Here you can find the metabolic models, implementation and results of the pipeline to infer community objectives by applying successive constraints and using OVA (Objective Variability Analysis).
- In the scritpt ```lacto_akker.m``` you can find the implementation for the objectives 0.5xBiomass(L) + 0.5xBiomass(A), 0.6xBiomass(A) + 0.4xBiomass(L) and for a flux point sampled from the Pareto front
- In the scritpt ```lacto_akker_Bio_ATPM.m``` you can find the implementation for the objective 0.5xBiomass(L) + 0.5xATPM(A).

The joint metabolic model is created and the fluxes sampled from the Pareto front are obtained in the ```generate_model_and_fluxes.ipynb``` notebook, which uses the MIMEco library. 

The resulting obejctive vectors can be found in the excel files, and the results of running OVA are stored in .mat objects, and can be loaded in the two scripts mentioned above. 
  

