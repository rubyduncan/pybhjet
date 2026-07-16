


def return_chains_eq_post_linear(path, summary=False): 

    '''returns a numpy array with each row corresponding to a step, with samples for each parameter at that step, 
    columns corresponding to each parameter 
    this is just the mean, median values that are being returned to be used 
    '''
    
    chains_path = Path(path + "/un_output")
    eq_wpost_path = chains_path / "chains/equal_weighted_post.txt" 
    equ_post_df = pd.read_csv(eq_wpost_path,sep='\s+',comment="#")

    equ_post_df.columns = [key_component_param(c) for c in equ_post_df.columns]
    labels = [latex_map.get(c, c) for c in equ_post_df.columns]

    if summary == True: 
        ndim = len(model_obj.free_parameters)  

        summary_eq = pd.DataFrame(
        {
            "name": equ_post_df.columns,
            "min": equ_post_df.min(axis=0).values,
            "max": equ_post_df.max(axis=0).values,
            "mean": equ_post_df.mean(axis=0).values,
            "median": equ_post_df.median(axis=0).values,
            "std": equ_post_df.std(axis=0).values,
        }
        )
        return summary_eq, chains_path
    
    else: 
        samples = equ_post_df.to_numpy(dtype=float)
        return samples, labels, chains_path

