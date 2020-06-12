import pystan
import pickle
import numpy

def check_div(fit):
    """Check transitions that ended with a divergence"""
    sampler_params = fit.get_sampler_params(inc_warmup=False)
    divergent = [x for y in sampler_params for x in y['divergent__']]
    n = sum(divergent)
    N = len(divergent)
    print('{} of {} iterations ended with a divergence ({}%)'.format(n, N,
            100 * n / N))
    if n > 0:
        print('Try running with larger adapt_delta to remove the divergences')

def check_treedepth(fit, max_depth = 10):
    """Check transitions that ended prematurely due to maximum tree depth limit"""
    sampler_params = fit.get_sampler_params(inc_warmup=False)
    depths = [x for y in sampler_params for x in y['treedepth__']]
    n = sum(1 for x in depths if x == max_depth)
    N = len(depths)
    print(('{} of {} iterations saturated the maximum tree depth of {}'
            + ' ({}%)').format(n, N, max_depth, 100 * n / N))
    if n > 0:
        print('Run again with max_depth set to a larger value to avoid saturation')

def check_energy(fit):
    """Checks the energy Bayesian fraction of missing information (E-BFMI)"""
    sampler_params = fit.get_sampler_params(inc_warmup=False)
    for chain_num, s in enumerate(sampler_params):
        energies = s['energy__']
        numer = sum((energies[i] - energies[i - 1])**2 for i in range(1, len(energies))) / len(energies)
        denom = numpy.var(energies)
        if numer / denom < 0.2:
            print('Chain {}: E-BFMI = {}'.format(chain_num, numer / denom))
            print('E-BFMI below 0.2 indicates you may need to reparameterize your model')

def _by_chain(unpermuted_extraction):
    num_chains = len(unpermuted_extraction[0])
    result = [[] for _ in range(num_chains)]
    for c in range(num_chains):
        for i in range(len(unpermuted_extraction)):
            result[c].append(unpermuted_extraction[i][c])
    return numpy.array(result)

def _shaped_ordered_params(fit):
    ef = fit.extract(permuted=False, inc_warmup=False) # flattened, unpermuted, by (iteration, chain)
    ef = _by_chain(ef)
    ef = ef.reshape(-1, len(ef[0][0]))
    ef = ef[:, 0:len(fit.flatnames)] # drop lp__
    shaped = {}
    idx = 0
    for dim, param_name in zip(fit.par_dims, fit.extract().keys()):
        length = int(numpy.prod(dim))
        shaped[param_name] = ef[:,idx:idx + length]
        shaped[param_name].reshape(*([-1] + dim))
        idx += length
    return shaped

def partition_div(fit):
    """ Returns parameter arrays separated into divergent and non-divergent transitions"""
    sampler_params = fit.get_sampler_params(inc_warmup=False)
    div = numpy.concatenate([x['divergent__'] for x in sampler_params]).astype('int')
    params = _shaped_ordered_params(fit)
    nondiv_params = dict((key, params[key][div == 0]) for key in params)
    div_params = dict((key, params[key][div == 1]) for key in params)
    return nondiv_params, div_params

def compile_model(filename, model_name=None, **kwargs):
    """This will automatically cache models - great if you're just running a
    script on the command line.

    See http://pystan.readthedocs.io/en/latest/avoiding_recompilation.html"""
    from hashlib import md5

    with open(filename) as f:
        model_code = f.read()
        code_hash = md5(model_code.encode('ascii')).hexdigest()
        if model_name is None:
            cache_fn = 'cached-model-{}.pkl'.format(code_hash)
        else:
            cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
        try:
            sm = pickle.load(open(cache_fn, 'rb'))
        except:
            sm = pystan.StanModel(model_code=model_code)
            with open(cache_fn, 'wb') as f:
                pickle.dump(sm, f)
        else:
            print("Using cached StanModel")
        return sm


def plot_rhat(filepath, filename):
    import pandas as pd
    import matplotlib.pyplot as plt
    print('filepath:', filepath)
    print('filename:', filename)

    mycsvfile=filepath+'/'+str(filename)+'.csv'
    print('mycsvfile:', mycsvfile)
    
    Summary =pd.read_csv(mycsvfile,  comment='#')
    #Summary =pd.read_csv(str(filename) + '.csv',  comment='#')

    hd_row=7                     
    Rhat=Summary['R_hat']
    Rhat=Rhat[hd_row:] 
    Rhat=pd.to_numeric(Rhat, errors='coerce')

    R_norm=sum(Rhat)/Rhat.size
    print('Normalized R=', float(R_norm))
    
    Rhat_large=((Rhat > 1.1).sum())
    print('Rhat >1.1: ', float(Rhat_large))

    plt.figure(figsize=(14,4))
    plt.subplot(1, 2, 1)
    Rhat.plot()
    plt.grid()
    plt.xlabel('Data points')
    plt.ylabel(r'$\hat R$')
    plt.subplot(1, 2, 2)
    pd.Series.hist(Rhat, bins=500)
    plt.text(1.005, 100.,  r'$R_{norm}:\ \ $'+str(R_norm) , size=10, color = 'red')
    plt.xlabel(r'$\hat R$ values')
    plt.ylabel('Count')
    plt.title(r'$\hat R$'); 
    plt.show()
    return


def plot_elbo(filename):
    with open(str(filename) + '.out' , 'r') as fd:
        counter=0
        lines=[]
        for line in fd.readlines():
            reading_line=line.strip().split(' ')
            if reading_line[0]=='Success!':
                start_iter=counter
            if reading_line[0]=='COMPLETED.':
                end_iter=counter
            counter+=1
            lines.append(line.strip().split(','))
        Convergence=[]  
        for row in arange(start_iter+4,(end_iter-5)+1):        
            line=lines[row][0]
            Convergence.append([x for x in re.split(',| ',line) if x!=''][0:2])
    Convergence=asarray(Convergence)  
    Iter=Convergence[:,0].astype('float')
    Elbo=Convergence[:,1].astype('float')
    
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    plt.plot(Elbo)
    my_yticks = ax.get_yticks()
    plt.yticks([])
    plt.xlabel('100*iter')
    plt.ylabel('ELBO')
    plt.show()
    return 