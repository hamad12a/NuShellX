      function radial(xn1,xl1,xn2,xl2,l)
cbab this is the radial me:   <xn1,xl1| r**l |xn2,xl2>
cbab in units of b**l
      implicit none
      integer,intent(in)::xn1,xn2,xl1,xl2,l
      real(kind=rc):: radial,xx,x2,x3xi,xj,x6p,x6,x7,x8,par,s,t,pi
      integer:: i,j,im,jm
      pi = 3.14159
      xx = l
      x2 = 2.0*xl1+2.0*xn1+1.0
      x3 = 2.0*xl2+2.0*xn2+1.0
      im = xn1+1.0
      jm = xn2 + 1.0
      s = 0.0
      do i = 1,im
        do j = 1,jm
           xi = i-1
           xj = j-1
           if(btest(xl1+xl2+xx),0) x6p = xi+xj+(xl1+xl2+xx+1.)/2.
           x6 = xl1+xl2+2.0*xi+2.0*xj+xx+1.0
           x7 = 2.0*xl1+2.0*xi+1.0
           x8 = 2.0*xl2+2.0*xj+1.0
           par=1.0
           if (btest(xi+xj,0) par=-par
           t = par/(factorial(xi)*factorial(xj))
           if(.not.btest(xl1+xl2+xx,0)) go to 31
30         t = t*factorial(x6p)
           xmm = (x6+1.)/2.
           go to 32
31         t = t*dfactorial(x6)
32         t = t/(dfactorial(x7)*dfactorial(x8)*factorial(xn1-xi)*factorial(xn2-xj))
           if(.not.btest(xl1+xl2+xx,0)) go to 13
40         t = t*(2.**xmm)/sqrt(pi)
13         s = s + t
        end do
      end do
      s = s*sqrt(dfactorial(x2)*dfactorial(x3))
      s = s*sqrt(factorial(xn1)*factorial(xn2))
      s = s/sqrt(2.0**(xn1+xn2+xx))
      radial = s
      return
      end
