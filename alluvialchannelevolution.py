def alluvialchannel(tend,uplift,kappa1,kappa2,deltaz):

    import numpy as np

    # Initial Topography
    nx=101
    dx=20
    xgrid=np.arange(0,nx*dx,dx)            # Grid
    area=np.zeros(nx)
    area=1000+0.5*xgrid**2                  # Hack's law relating drainage area and stream length

    # Define Topography Based on Slope
    #topo=np.zeros(nx)
    #for i in range(nx-2,-1,-1):
    #    topo[i]=topo[i+1]+0.0175*dx
        
    #topo=topo-min(topo);
    
    topo=np.loadtxt('steadystatealluvial_U0pt02.txt')

    slope=np.zeros(nx)
    kappa=np.ones(nx)
    halfway=round(0.5*nx)
    kappa[0:halfway]=kappa1*kappa[0:halfway]
    kappa[halfway:nx]=kappa2*kappa[halfway:nx]
    topoold=np.zeros(nx)
    flux=np.zeros(nx+1)
    
    area0=area[0]
    slope0=0.1

    t=0
    dt=0.0005
    while t<tend:
    
        topoold=topo[:]
    
        slope[0:nx-1]=1/dx*abs((topo[1:nx]-topo[0:nx-1]))
        slope[nx-1]=slope[nx-2]
        
        flux[1:nx+1]=kappa*(slope*area**(4/10))**(3/2)
        flux[0]=kappa[0]*(slope0*area0**(4/10))**(3/2)
        
        effectivemaxD=max((kappa*slope*area**(4/10))**(3/2)/slope)
    
        dt=0.1*dx*dx/effectivemaxD
    
        topo=topo+dt*uplift-dt/dx*(flux[1:nx+1]-flux[0:nx])
        topo[nx-1]=topoold[nx-1]
    
        if (t>0.0001) and (deltaz>0):
            topo[nx-1]=topo[nx-1]-deltaz
            deltaz=0;
    
        t=t+dt

    return(xgrid,area,topo,slope)


