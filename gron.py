#
import numpy as np
import matplotlib.pyplot as plt
#
#   objective
#
def obj(x):
#
    a = 1.
    b = 32.
    c = -1.
    f = a*np.sin(b*x)*np.exp(c*x)
    df = a*c*np.sin(b*x)*np.exp(c*x) + a*np.exp(c*x)*np.cos(b*x)*b
    ddf = a*c*b*np.cos(b*x)*np.exp(c*x) + a*c*c*np.sin(b*x)*np.exp(c*x)  + a*c*np.exp(c*x)*np.cos(b*x)*b - a*np.exp(c*x)*np.sin(b*x)*b*b
#
    return [f, df, ddf]
#
#   constraint
#
def con(x):
#
    a = 1./4.
    b = 16
#
    g = a*np.cos(b*x)
    dg = -a*b*np.sin(b*x)
    ddg = -a*b*b*np.cos(b*x)
#
    return [g, dg, ddg]
#
#   main
#
if __name__ == '__main__':
#
#   prelims
#
    x_dom = np.linspace(0.4,0.6,num=1000)
    x_app_dom = np.linspace(0,1,num=1000)
    f_dom = np.zeros_like(x_dom) 
    g_dom = np.zeros_like(x_dom) 
#
    for i in range(len(x_dom)):
        f_dom[i] = obj(x_dom[i])[0]
        g_dom[i] = con(x_dom[i])[0]
#
    x = 0.5
#
    cg = 1; cf = 1
    f_app=0.; g_app=0.
    mov=0.1
    xnew=0.
    fold=1e8
    xold=0
#
#   loop
#
    for k in range(100):
#
#       simulate       
#
        [f,df,ddf] = obj(x)
        [g,dg,ddg] = con(x)
#
#       check if feasible descent, otherwise check if approximations are conservative
#
        flg=0
        if f > fold or g > 0:
            if k>0 and f_app < f or g_app < g:
                with open('tmp_%d.tex'%(k-1),'a') as file:
                    file.write('-> not accepted\n\n')
                x=xold
                if f_app < f:
                    cf = cf*2.
                    with open('tmp_%d.tex'%(k-1),'a') as file:
                        file.write('--> objective\n\n')
                if g_app < g:
                    cg = cg*2.
                    with open('tmp_%d.tex'%(k-1),'a') as file:
                        file.write('--> constraint\n')
                [f,df,ddf] = obj(x)
                [g,dg,ddg] = con(x)
            else:
                with open('tmp_%d.tex'%(k-1),'a') as file:
                    file.write('-> accepted\n\n')
                flg=1
                cg = 1
                cf = 1
        else:
            with open('tmp_%d.tex'%(k-1),'a') as file:
                file.write('-> feasible descent\n\n')
            flg=1
#
        print(k,flg, abs(x-xold))
#
        if flg == 1 and abs(x-xold)<1e-6:
            exit()
#
        fold=f
#
#       enforce strict convexity
#
        ddf=max(ddf,1e-6)
        ddg=max(ddg,0.)
#
#       plot the approximations around the current point
#
        f_app_dom = f + (x_app_dom-x)*df + cf*0.5*ddf*(x_app_dom-x)**2.
        g_app_dom = g + (x_app_dom-x)*dg + cg*0.5*ddg*(x_app_dom-x)**2.
#
#       QPQC update
#
        lab_lo = 1e-9
        lab_up = 1e9
        while (lab_up-lab_lo)/(lab_lo+lab_up)>1e-9:
            lab=0.5*(lab_up+lab_lo)
            xdel=np.minimum(np.maximum((df+lab*dg)/(ddf*cf+lab*ddg*cg),-mov),mov)
            xnew=np.minimum(np.maximum(0.,x-xdel),1.)
            gt=g+dg*(xnew-x)+cg*ddg*(xnew-x)*(xnew-x)/2.
            if gt>0 :
                lab_lo=lab
            else:
                lab_up=lab
#
#       get the approximate function values at new point
#
        f_app = f + (xnew-x)*df + cf*0.5*ddf*(xnew-x)**2.
        g_app = g + (xnew-x)*dg + cg*0.5*ddg*(xnew-x)**2.
        [f_new,_,_] = obj(xnew)
        [g_new,_,_] = con(xnew)
#
#       plot various
#
        plt.plot(x_dom,f_dom,'k--' ,x_dom,g_dom,'k:', \
            x_app_dom,f_app_dom,'b-', x_app_dom,g_app_dom,'r-',\
            [x],[f],'bo',[x],[g],'ro', \
            [xnew], [g_app], 'r*', [xnew], [f_app], 'b*', \
            [xnew], [g_new], 'rx', [xnew], [f_new], 'bx', [x], [0], 'k|')
        plt.grid()
        ax = plt.gca()
        ax.set_xlim([0.4, 0.6])
        ax.set_ylim([-1, 1])
        txt = 'yoyo'
#
        plt.tight_layout()
#
        with open('tmp_%d.tex'%k,'w') as file:
            txt='-'*30
            file.write(txt+'\n\n')
            txt='Iteration %d'%k
            file.write(txt+'\n\n')
            txt='-'*30
            file.write(txt+'\n\n')
            txt='l = %f'%lab
            file.write(txt+'\n\n')
            txt= 'x = %f'%x
            file.write(txt+'\n\n')
            txt= 'f = %f'%f
            file.write(txt+'\n\n')
            txt= 'df = %f'%df
            file.write(txt+'\n\n')
            txt= 'ddf = %f'%ddf
            file.write(txt+'\n\n')
            txt= 'cf = %f'%cf
            file.write(txt+'\n\n')
            txt= 'g = %f'%g
            file.write(txt+'\n\n')
            txt= 'dg = %f'%dg
            file.write(txt+'\n\n')
            txt= 'ddg = %f'%ddg
            file.write(txt+'\n\n')
            txt= 'cg = %f'%cg
            file.write(txt+'\n\n')
            txt= 'xnew = %f'%xnew
            file.write(txt+'\n\n')
            txt= 'fa[xnew] = %f'%f_app
            file.write(txt+'\n\n')
            txt= 'ga[xnew] = %f'%g_app
            file.write(txt+'\n\n')
            txt= 'f[xnew] = %f'%f_new
            file.write(txt+'\n\n')
            txt= 'g[xnew] = %f'%g_new
            file.write(txt+'\n\n')
#
        plt.savefig('itr_%d.png'%k)
        plt.close()
#
#       update x
#
        xold=x
        x=xnew
#

