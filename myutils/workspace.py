import os
import sys
from pathlib import Path

def workspace(mydir:str='Selina',project:str='project',script:str | None=None) -> tuple:
    """get workspace paths for script

    Args:
        mydir (str, optional): Name of personal directory (not a path). Defaults to 'Selina'.
        project (str, optional): Name of project directory in mydir (not a path). Defaults to 'project'.
        script (str, optional): Give a path to the script file. Defaults to None.

    Returns:
        tuple: (path_to_workspace,path_to_experiment)
    """
    
    if script!=None:
        workspace = next(p for p in Path(script).resolve().parents if p.name == mydir)/project
        experiment = next(p for p in Path(script).resolve().parents if "experiment" in p.name.lower())
    else: 
        cwd = Path.cwd()
        workspace = next(p for p in cwd.resolve().parents if p.name == mydir)/project
        experiment = workspace

        is_notebook = False
        from IPython import get_ipython
        if get_ipython() is not None: is_notebook = True

        if not is_notebook:
            cwd =  Path(sys._getframe(1).f_code.co_filename).parent
        if any("experiment" in p.name.lower() for p in cwd.resolve().parents):
            experiment = next(p for p in cwd.resolve().parents if "experiment" in p.name.lower())

    return workspace, experiment
