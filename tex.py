def tex(k,lab,x,f,df,ddf,cf,g,dg,ddg,cg,fapp,gapp,fnew,gnew,xnew):
        with open('tmp2_%d.tex'%k,'w') as file:
            txt='QP subproblem at $\\x^{k}$:'
            file.write(txt+'\n\n')
            file.write('\\bigskip\n')
            txt= '$\\x^k= %7.3f$'%x
            file.write(txt+'\n\n\n')
            txt= '$\\f^k = %7.3f$'%f
            file.write(txt+'\n\n')
            txt= '$\\d^k_{\\x} \\f = %7.3f$'%df
            file.write(txt+'\n\n')
            txt= '$\\d\\d^k_{\\x} \\tilde \\f = %7.3f$'%ddf
            file.write(txt+'\n\n')
            txt= '$\\alpha_{\\footnotesize\\f} = %7.3f$'%cf
            file.write(txt+'\n\n')
            file.write('\\bigskip\n')
            txt= '$\\g^k = %7.3f$'%g
            file.write(txt+'\n\n')
            txt= '$\\d^k_{\\x}\\g = %7.3f$'%dg
            file.write(txt+'\n\n')
            txt= '$\\d\\d^k_{\\x} \\tilde \\g = %7.3f$'%ddg
            file.write(txt+'\n\n')
            txt= '$\\alpha_{\\footnotesize\\g} = %7.3f$'%cg
            file.write(txt+'\n\n')
            file.write('\\bigskip\n')
            txt='with Lagrange multiplier'
            file.write(txt+'\n\n')
            txt='$\\lambda^{k\\star} = %7.3f$'%lab
            file.write(txt+'\n\n')
            txt='at the solution'
            file.write(txt+'\n\n')
            txt= '$\\x^{k\\star} = %7.3f$'%xnew
            file.write(txt+'\n\n')
        with open('tmp1_%d.tex'%k,'w') as file:
            txt='At the solution of the QP the function approximations have the values'
            file.write(txt+'\n\n')
            txt= '$\\tilde \\f^{k\\star} = %13.6e$'%fapp
            file.write(txt+'\n\n')
            txt= '$\\tilde \\g^{k\\star} = %13.6e$'%gapp
            file.write(txt+'\n\n')
            file.write('\\bigskip\n')
            txt='while the actual functions are evaluated to be'
            file.write(txt+'\n\n')
            txt= '$\\f^{k\\star} = %13.6e$'%fnew
            file.write(txt+'\n\n')
            txt= '$\\g^{k\\star} = %13.6e$'%gnew
            file.write(txt+'\n\n')
