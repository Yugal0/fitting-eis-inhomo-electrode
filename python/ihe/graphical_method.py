
def generate_graph(model_name,R1_range,R2_range):

    """generate contour plots of theta_dip with tau1 and tau2 for the graphical method of fitting[1]
       Choudhary et al. 2023, J-ECS, Page number x to y
    
    Parameters
    ----------
    model_name : str
        name of the model, can be "TLMtwo", "TLMlin", "TLMilin" 
    R1_range : float
        limit on value of R1 in the plot
    R2_range : float
        limit on value of R2 in the plot
    """