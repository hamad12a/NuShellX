      function rl(i1,i2,l)
cbab this is the radial me:   <n(i1)| r**l |n(i2)>
cbab in units of b**l
      common /indexm/ xl(3,200)
      dfacn(x) = dfac(x)
      pi = 3.14159
      xn1 = xl(1,i1)-1.0
      xl1 = xl(2,i1)
      xn2 = xl(1,i2)-1.0
      xl2 = xl(2,i2)
      xx = l
      x2 = 2.0*xl1+2.0*xn1+1.0
      x3 = 2.0*xl2+2.0*xn2+1.0
      im = xn1+1.0
      jm = xn2 + 1.0
      s = 0.0
      do 13 i = 1,im
      do 13 j = 1,jm
      xi = i-1
      xj = j-1
      if(par(xl1+xl2+xx)) 20,20,21
20    x6p = xi+xj+(xl1+xl2+xx+1.)/2.
21    continue
      x6 = xl1+xl2+2.0*xi+2.0*xj+xx+1.0
      x7 = 2.0*xl1+2.0*xi+1.0
      x8 = 2.0*xl2+2.0*xj+1.0
      t = par(xi+xj)/(fac(xi)*fac(xj))
      if(par(xl1+xl2+xx)) 30,30,31
30    t = t*fac(x6p)
      xmm = (x6+1.)/2.
      go to 32
31    t = t*dfacn(x6)
32    continue
      t = t/(dfacn(x7)*dfacn(x8)*fac(xn1-xi)*fac(xn2-xj))
      if(par(xl1+xl2+xx)) 40,40,13
40    t = t*(2.**xmm)/sqrt(pi)
13    s = s + t
      s = s*sqrt(dfacn(x2)*dfacn(x3))
      s = s*sqrt(fac(xn1)*fac(xn2))
      s = s/sqrt(2.0**(xn1+xn2+xx))
      rl = s
      return
      end
