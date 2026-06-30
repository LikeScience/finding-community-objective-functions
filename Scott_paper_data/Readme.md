# Recovering objectives in _C. autoethanogenum & C. kluyveri community_

Here you can find the metabolic models, implementation and results of running invFBA to to infer community objectives, generating fluxes through forward FBA. 
- In the scritpt ```auto_kluy_invFBA.m``` you can find the implementation for the objectives 0.5xBiomass(A) + 0.5xBiomass(K)
- In the script ```auto_kluy_invFBA_OVA.m``` inside the subdirectory ```OVA``` you can find the implementation of a pipeline using OVA to reduce solution space and the results of running this pipeline  for the objective 0.5xBiomass(L) + 0.5xATPM(A).

The joint metabolic model is created and the fluxes sampled from the Pareto front are obtained in the ```generate_model_and_fluxes.ipynb``` notebook, which uses the MIMEco library. 

The resulting obejctive vectors can be found in the excel files, and the results of running OVA are stored in .mat objects, and can be loaded in the two scripts mentioned above. 
  

