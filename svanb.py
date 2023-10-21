#
import logging
import numpy as np
import matplotlib.pyplot as plt
logging.getLogger('matplotlib').setLevel(level=logging.CRITICAL)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica",
    "font.serif": "Times New Roman"
})
#
from tex import tex
from fun import obj,con
#
#   main
#
if __name__ == '__main__':
#
#   prelims and settings
#
    x_dom = np.linspace(0.4,0.6,num=1000)
    xapp_dom = np.linspace(0,1,num=1000)
    f_dom = np.zeros_like(x_dom) 
    g_dom = np.zeros_like(x_dom) 
    fapp=0.; gapp=0.
    cg = 1; cf = 1
    xnew=0.0
    xold=0.0
    mov=0.1
    eps=1e-6
#
#   for plotting
#
    for i in range(len(x_dom)):
        f_dom[i] = obj(x_dom[i])[0]
        g_dom[i] = con(x_dom[i])[0]
#
#   starting point
#
    x = 0.5
#
    print('%3s %14s %14s %14s %14s %14s'%('k','f','g','x','lab','abs(x-xold)'))
#
#   loop
#
    for k in range(99):
#
#       simulate / evaluate functions       
#
        [f,df,ddf] = obj(x)
        [g,dg,ddg] = con(x)
#
#       check if approximations are conservative
#
        flg=0
        if k>0 and fapp < f-eps or gapp < g-eps:
            with open('tmp1_%d.tex'%(k-1),'a') as file:
                file.write('\\bigskip\n However, the solution (step) \
                    is deemed \\emph{unacceptable}, because \n \n')
            x=xold
            if fapp < f-eps:
                cf = cf*2.
                with open('tmp1_%d.tex'%(k-1),'a') as file:
                    file.write('$\\to$ the value of the objective \
                        function approximation is less than the actual function \
                        value, at the new design point, \
                        $\\tilde {\\f}^{k\\star} < {\\f}^{k\\star}$.')
            if gapp < g-eps:
                cg = cg*2.
                with open('tmp1_%d.tex'%(k-1),'a') as file:
                    file.write('\n\n $\\to$ the value of the constraint function \
                        approximation is less than the actual function value, at \
                        the new design point, $\\tilde {\\g}^{k\\star} < {\\g}^{k\\star}$.')
            with open('tmp1_%d.tex'%(k-1),'a') as file:
                file.write('\n\n \\bigskip \n\n That is, at least one approximation \
                    function is not \\emph{conservative}.\n')
            [f,df,ddf] = obj(x)
            [g,dg,ddg] = con(x)
        else:
            with open('tmp1_%d.tex'%(k-1),'a') as file:
                file.write('\\bigskip \n Both the value objective function approximation \
                    and the constraint function approximation, at the new design point \
                    $x^{k\\star}$, is \\emph{conservative} with respect to the actual function \
                    values, $\\tilde {\\f}^{k\\star} > {\\f}^{k\\star}$, \
                    $\\tilde {\\g}^{k\\star} > {\\g}^{k\\star}$ \n\n')
            flg=1; cg=1; cf=1
#
        if flg == 1 and abs(x-xold)<eps:
            with open('tmp1_%d.tex'%(k-1),'a') as file:
                file.write('\n\n\\bigskip Terminated on $|\\x^{k\\star}-\\x^{k}|$$<$%8.1e\n'%eps)
            exit()
#
#       plot the approximations around the current point
#
        fapp_dom = f + (xapp_dom-x)*df + cf*0.5*ddf*(xapp_dom-x)**2.
        gapp_dom = g + (xapp_dom-x)*dg + cg*0.5*ddg*(xapp_dom-x)**2.
#
#       enforce strict convexity
#
        ddf=max(ddf,eps)
        ddg=max(ddg,eps)
#
#       QPQC update
#
        lab_lo = 1e-9
        lab_up = 1e9
        while (lab_up-lab_lo)/(lab_lo+lab_up)>eps:
            lab=0.5*(lab_up+lab_lo)
            xdel=min(max((df+lab*dg)/(ddf*cf+lab*ddg*cg),-mov),mov)
            xnew=min(max(0.,x-xdel),1.)
            gt=g+dg*(xnew-x)+cg*ddg*(xnew-x)*(xnew-x)/2.
            if gt>0 :
                lab_lo=lab
            else:
                lab_up=lab
#
#       get the approximate function values at new point
#
        fapp = f + (xnew-x)*df + cf*0.5*ddf*(xnew-x)**2.
        gapp = g + (xnew-x)*dg + cg*0.5*ddg*(xnew-x)**2.
        [fnew,_,_] = obj(xnew)
        [gnew,_,_] = con(xnew)
#
#       plot various
#
        plt.plot(x_dom,f_dom,'k--',label='$\\textrm{F}$') 
        plt.plot(x_dom,g_dom,'k:',label='$\\textrm{G}$') 
        plt.plot(xapp_dom,fapp_dom,'b-',label='$\\tilde \\textrm{F}$') 
        plt.plot(xapp_dom,gapp_dom,'r-',label='$\\tilde \\textrm{G}$') 
        plt.plot([x],[f],'bo',label='$\\tilde \\textrm{F}^k$') 
        plt.plot([x],[g],'ro',label='$\\tilde \\textrm{G}^k$') 
        plt.plot([xnew],[fapp],'b*',label='$\\tilde \\textrm{F}^{k\\star}$') 
        plt.plot([xnew],[fnew],'bx',label='$\\textrm{F}^{k\\star}$') 
        plt.plot([xnew],[gapp],'r*',label='$\\tilde \\textrm{G}$') 
        plt.plot([xnew],[gnew],'rx',label='$\\textrm{G}^{k\\star}$') 
        plt.plot([x],[0],'k|') 
        plt.legend(loc="lower left",ncol=5)
        plt.xlabel("$x$")
        plt.grid()
        ax = plt.gca(); ax.set_xlim([0.45, 0.6]); ax.set_ylim([-1.5, 1])
        plt.tight_layout()
        plt.savefig('itr_%d.eps'%k,format='eps')
        plt.close()
#
#       screen output
#
        print('%3d %14.3e %14.3e %14.3e %14.3e %14.3e'%(k,f,g,x,lab,abs(x-xold)))
#
#       tex output
#
        tex(k,lab,x,f,df,ddf,cf,g,dg,ddg,cg,fapp,gapp,fnew,gnew,xnew)
#
#       update x
#
        xold=x
        x=xnew
#

