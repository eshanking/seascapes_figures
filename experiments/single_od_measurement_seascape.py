from fears.utils import AutoRate

folder_path = '/Users/kinge2/repos/fears/fears/data/od_plates_no_lid'

layout_path = '/Users/kinge2/repos/fears/fears/data/plate_layout.csv'
reference_path = '/Users/kinge2/repos/fears/fears/data/plates/20210929_plate1.csv'

t_obs = 12*3600

e = AutoRate.Experiment(folder_path=folder_path,mode='single_measurement',
                        exp_layout_path=layout_path,
                        ref_data_path=reference_path,t_obs=t_obs)
e.execute()