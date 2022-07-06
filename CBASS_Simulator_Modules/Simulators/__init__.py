from CBASS_Simulator_Modules.Simulators import Ground, Noise, Sky, Load,Dipole

model_functions = {'ground':Ground.GroundModel,
                   'dipole':Dipole.DipoleModel,
                   'fnoise':Noise.FnoiseModel,
                   'wnoise':Noise.WnoiseModel,
                   'sky'   :Sky.SkyModel,
                   'load'  :Load.LoadModel}
