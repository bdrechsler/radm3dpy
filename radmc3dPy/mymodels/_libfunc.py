import subprocess as sp

def updateModelList():
    """Updates the list of models in the library
    
    """
    dum = sp.Popen(['ls -1 '+mod_path+'/models/*.py'], shell=True, stdout=sp.PIPE, stderr=sp.PIPE).communicate()[0].split()





def getModelNames():
    """Returns the name of the available models

    """

    mod_names = []
   
    # Get the name of all model files in the module directory
    import radmc3dPy
    mod_path = radmc3dPy.__file__.strip()[:-12]
    dum = sp.Popen(['ls -1 '+mod_path+'/models/*.py'], shell=True, stdout=sp.PIPE, stderr=sp.PIPE).communicate()[0].split()

    for i in range(len(dum)):
        sdum = dum[i].split('/')
        #id1 = dum[i].index('model_')+6
        #id1 = dum[i].index('.py')
        #modname = dum[i][id1:id2]
        modname = sdum[len(sdum)-1][:-3]

        if ((modname!='template')&(modname!='__init__')&(modname!='_libfunc')):
            mod_names.append(modname)

    ## Get the name of all model files in the current working directory
    #if os.getcwd().strip()!=mod_path:
        #mod_path = os.getcwd()
        #dum = sp.Popen(['ls -1 '+mod_path+'/model_*.py'], shell=True, stdout=sp.PIPE, stderr=sp.PIPE).communicate()[0].split()

        #for i in range(len(dum)):
            #id1 = dum[i].index('model_')+6
            #id2 = dum[i].index('.py')
            #modname = dum[i][id1:id2]
            #if modname!='template':
                #mod_names.append(modname)

    return mod_names
    

def getModelDesc(model=''):
    """Returns the brief description of the selected model

    """

    if model.strip()=='':
        print 'ERROR'
        print 'No model name is given'
        return

    try:
        mdl = __import__(model)
    except:
        try:
            mdl  = __import__('radmc3dPy.models.'+model, fromlist=['']) 
        except:
            print 'ERROR'
            print ' model_'+model+'.py could not be imported'
            print ' The model files should either be in the current working directory or'
            print ' in the radmc3d python module directory'
            return -1

    
    if callable(getattr(mdl, 'getModelDesc')):
        return mdl.getModelDesc()
    else:
        print 'ERROR'
        print 'model_'+model+'.py does not contain a getModelDesc() function.'
        return






