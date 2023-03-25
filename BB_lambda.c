#include <stdio.h>
#include <math.h>
#include <cpgplot.h>
#include <string.h>
#include <stdlib.h>
 
// gcc -c BB_lambda.c -I/usr/local/pgplot -o BB_lambda.o; gfortran -o BB_lambda BB_lambda.o -I/opt/local/include -L/opt/local/lib -lcpgplot -lpgplot -lpng -lz; ./BB_lambda

int main()
{
  int i, censor, q, numChan, n = 1000; 
  float  x1,x2,x3,y1,y2,y3,xmin, xmax, ymin, ymax, min, peak_nu; 
  float  xline[n], yline[n],zline[n], m[n], T, R;
  FILE *f1,*f2;
  char  format[10],  star[10], lineInput[1000], systemCall[100], ps_out[100], ps_out_dev[100];
 
  float ybottom = -8, ytop = 17; 
   float xbottom = -8, xtop = 0.2;//log

 double I, B, h = 6.62606957e-34, K = 1.3806488e-23, L, power[n], *pee, P, intensity[n],  rate[n], *ree, RATE; 
  float lambda_1, lambda_2, c = 299792458, pi =  3.14159;
  
  sprintf(ps_out_dev, "BB_lambda.eps/cps");
  sprintf(ps_out, "BB_lambda.eps");

  printf("Doing for Sun [s] or other star [o]? "); scanf("%s", star); 
  //strcpy(star,"s");
  if (strcmp(star, "s") == 0) {T = 5780; R = 696300e3;}
  else{
    printf("Surface temperature of star in K? ");   scanf("%f", &T);
    printf("Radius of star in metres? ");   scanf("%f", &R); //T = 3e4, R = 5e9; //Zdenka p103
  } // SEE 343_Applied_Physics/Assignments+Tests/ BB_star-freq.c  FOR OTHER STAR
  
  
  printf("Output format, screen [x] or file [p]?");
  scanf("%s", format);
  if (strcmp(format, "x") == 0 ) cpgopen("/xs");
  else  cpgopen(ps_out_dev);
  cpgpap(8.0,0.7); 
   
  cpgsch(1.6);
  cpgqvp(0,&x1,&x2,&y1,&y2);
  cpgsvp(x1+0.07,x2,y1+0.05,y2-0.05);
  cpgswin(xbottom, xtop, ybottom, ytop); 
 
  cpgslw(3);  
 
 cpglab("","Brightness, log\\d10\\uB\\d\\gl\\u [W m\\u-2\\d m\\u-1\\d sr\\u-1\\d]","");
  cpgmtxt("B",2.5,0.5,0.5, "Wavelength, \\gl [m]");
  cpgmtxt("T",2.2,0.5,0.5, "Frequency, \\gn [Hz]");
  //////////////////////////////////////////////////////////      
 
  int style = 0; //0;
  int colour = 1;
  cpgsch(1.3);
  cpgtext(-4, 14.3,"Planck function for T = 5780 K");
  cpgsch(1.6);
  

 cpgslw(10);  // ADD RAINBOW
    float lam = 380e-9;
    float bri = 2*h*pow(c,2)/(pow(lam,5)*(exp((h*c)/(lam*K*T))-1));
    cpgsci(12);cpgmove(log10(lam), ybottom); cpgdraw(log10(lam), log10(bri)); //violet

    lam = 450e-9;bri = 2*h*pow(c,2)/(pow(lam,5)*(exp((h*c)/(lam*K*T))-1));
    cpgsci(4);cpgmove(log10(lam), ybottom); cpgdraw(log10(lam), log10(bri)); //blue

    lam = 520e-9;bri = 2*h*pow(c,2)/(pow(lam,5)*(exp((h*c)/(lam*K*T))-1));
    cpgsci(3);cpgmove(log10(lam), ybottom); cpgdraw(log10(lam), log10(bri)); //green

    lam = 590e-9;bri = 2*h*pow(c,2)/(pow(lam,5)*(exp((h*c)/(lam*K*T))-1));
    cpgsci(7);cpgmove(log10(lam), ybottom); cpgdraw(log10(lam), log10(bri)); //yellow

    lam = 660e-9;bri = 2*h*pow(c,2)/(pow(lam,5)*(exp((h*c)/(lam*K*T))-1));
    cpgsci(8);cpgmove(log10(lam), ybottom); cpgdraw(log10(lam), log10(bri)); //orange
    
    lam = 730e-9;bri = 2*h*pow(c,2)/(pow(lam,5)*(exp((h*c)/(lam*K*T))-1));
    cpgsci(2);cpgmove(log10(lam), ybottom); cpgdraw(log10(lam), log10(bri)); //red

    cpgslw(3); cpgsci(1);
  
   
 //////////////////  TICKS, START WITH LAMBDA ALONG BOTTOM ////////
  float tick[1000], subtick[10],main_tick_size = 0.5, sub_tick_size = 0.25, sub_space = 5, tickstart = log10(3e8/pow(xtop,10)), tickend = log10(3e8/pow(xbottom,10)); // converting to freq
  char step_string[20];
  float tick_label_off = -0.03;
  int step = 1000;
 
  int tick_start, tick_end;
  
 
 ///////////////// HATCHING ////////////////////////////////////////////////////
  float RJ = (h*c)/(K*T);  // want to shade this
  cpgsci(15);cpgsfs(1); //cpgshs(45,3,0); //changing hatching - angle, space, something
  //cpgrect(log10(RJ),xtop,ybottom,ytop); 
  ////////////////////////////////////////////////////////////////////////////

 cpgsci(1);

   int j;
   float nu_start = log(c) - xtop, nu_end = log(c) - xbottom;
   for (q=xbottom ; q<=xtop ; q++){
     tick[q] = log10(c) - q ;
     float pos_m = (xtop-q)/(xtop-xbottom); 
     cpgtick(xbottom, ytop, xtop, ytop, pos_m, 0.5, 0.5, 0, 0, ""); 
     sprintf(step_string,"10\\u%1.0f",tick[q]);tick_label_off = 0;
     cpgtext(q - tick_label_off, ytop + (ytop-ybottom)/40,  step_string);
     // try subticks here 
     for (j = 0; j <= 9; j++){
       float pos_s =   (xtop - (q+ log10(1+j)))/(xtop-xbottom);// - subtick[j]; 
       cpgtick(xbottom, ytop, xtop, ytop, pos_s, 0.2, 0.2, 0, 0, "");
     }
   }

  for (i = 0; i < n; i++) {
    xline[i] =  xbottom + ((xtop - xbottom)/n)*i;  
    zline[i] = pow(10,xline[i]);  
    B = 2*h*pow(c,2)/((pow(zline[i],5)*(exp((h*c)/(zline[i]*K*T))-1))); // in lambda
    //I = pi*B;
    yline[i] = log10(B); 
   /////shading of ionis part
    cpgsci(15);cpgslw(3);
    /* if (xline[i-1] < log10(912e-10)){ */
/*       cpgmove(xline[i-1],ybottom);  */
/*       cpgdraw(xline[i-1],yline[i-1]); */
/*     } */
    // printf("%1.2f to  %1.2f, wave= %1.2e m => %1.2e Hz, I = %1.2e  yline[%d] =%1.2f\n", xbottom, xtop, zline[i],  c/zline[i], I, i, yline[i]);
  }
     cpgsls(1); cpgsci(1);cpgslw(7);
    cpgline(n,xline,yline);
   
cpgslw(3); cpgsci(1);
cpgbox("bcnl",0.,0.,"bcsnt",0.,0.); // over hatching
 for (q=xbottom ; q<=xtop ; q++) cpgtick(xbottom, ybottom, xtop, ybottom, (q-xbottom)/(xtop-xbottom), 0.5, 0.5, 0, 0, ""); // main ticks
   // subticks
   for (tick_start = xbottom-1; tick_start <= xtop+1;  tick_start++){ // have to loop from -3 to -2, -2 to -1,
     tick_end = tick_start+ 1;
     for (q=0 ; q<=10; q++) { // 10 per main ticksubticks 
       tick[q] = tick_start +  log10(q)/(tick_end - tick_start);
       cpgtick(tick_start, ybottom, tick_end, ybottom, (tick[q]-tick_start)/(tick_end - tick_start), 0.2, 0.2, 0, 0, "");
     } 
   }
   // now down the way for freq
  
    if (strcmp(format, "x") != 0) sprintf(systemCall,"gv %s &", ps_out); system(systemCall);
  
   cpgend();

}
