// Jordan Ball
// November 2016
// Reactor Physics Project
// Coding Language C++

// This script can be used along with a neutron beam and detector in order to locate sub-surface oil repositories.

// Everything kept in cm

#include <iostream>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <math.h>

using namespace std;

int main()
{
    double Eo = 2500000;
    double E;
    double d = 1;
    double z = 0;
    double No;
		
    double ro_stone = 2.2+0.01*z/1000;          // in g/cm^3
    double ro_oil = 0.975;                      // in g/cm^3
    double ro;
	
    double Na = 6.022*pow(10,23);
    double A_stone = 100;
    double A_oil = 28.231;                      // nonphysical approximation
    double alpha_stone = 99/100;
    double alpha_oil = 27.231/29.231;
    double A;
    double alpha;
	
    srand(740982364);                           //((double)rand()/(double)RAND_MAX) generates random number between 0 and 1
	
    double stone_slow_scat = 19.65361*pow(10,-24);
    double stone_slow_abs = 0.168101*pow(10,-24);
    double stone_fast_scat = 15.95968*pow(10,-24);
    double stone_fast_abs = 3.821623*pow(10,-24);
    double oil_slow_scat = 83.84325*pow(10,-24);
    double oil_slow_abs = 0.010126*pow(10,-24);
    double oil_fast_scat = 36.49341*pow(10,-24);
    double oil_fast_abs = 22.61355*pow(10,-24);
    double abs;
    double scat;
    double sum_macro;	
	
    double r;
    int i;
    double interact = 0;
    double absorb = 0;
	
    E = Eo;

    cout<<"Input initial neutron number:"<<endl;
    cin>>No;
	
    int N = No;
    int reflected = 0;
	
    for (i = 0; i < No; i++){       //1 cycle for each neutron
		
        z = 1/(ro_stone*(stone_slow_abs+stone_slow_scat)*Na/A_stone);       //starting it with one MFP so the code doesn't end up with one less collision per neutron	
			
        do{			//run this until absorbed or past 40m
			
            if (z <= 0){        //removes neutrons above the surface as prob. interaction is negligible in comparison
                reflected = reflected + 1;
                N = N-1;
                z=50;
            }
		
            else {      //For all neutrons not reflected
			
                ro_stone = 2.2+0.01*z/1000;
			
			
				
                if (z <= 20){       //Splits the neutrons by medium
					
                    ro = ro_stone;
                    alpha = alpha_stone;
                    A = A_stone;
				
                    if (E <= 1){    //Splits by energy group
					
                        abs = stone_slow_abs;
                        scat = stone_slow_scat;				
                    }
					
                    else {
					
                        abs = stone_fast_abs;
                        scat = stone_fast_scat;
                    }
                }
			
                else if (20 < z <= 40) {
				
                    if (E < 1){
					
                        r = ((double)rand()/(double)RAND_MAX);
                        sum_macro = 0.3*ro_oil*(oil_slow_abs+oil_slow_scat)*(Na/A_oil) + 0.7*ro_stone*(stone_slow_abs+stone_slow_scat)*(Na/A_stone);
                        
                            if (r <= 0.3*ro_oil*(oil_slow_abs+oil_slow_scat)*(Na/A_oil)/sum_macro){		//To determine if it interacts with oil or stone
                            
                                abs = oil_slow_abs;
                                scat = oil_slow_scat;
                                ro = ro_oil;
                                alpha = alpha_oil;
                            A = A_oil;
                        }
						
                        else {
						
                            abs = stone_slow_abs;
                            scat = stone_slow_scat;
                            ro = ro_stone;
                            alpha = alpha_stone;
                            A = A_stone;
                        }	
                    }
				
                    else {
					
                        r = ((double)rand()/(double)RAND_MAX);
                        sum_macro = 0.3*ro_oil*(oil_fast_abs+oil_fast_scat)*Na/A_oil + 0.7*ro_stone*(stone_fast_abs+stone_fast_scat)*Na/A_stone;
                                            
                        if (r <= (0.3*ro_oil*(oil_fast_abs+oil_fast_scat)*Na/A_oil)/sum_macro){
						
                            abs = oil_fast_abs;
                            scat = oil_fast_scat;
                            ro = ro_oil;
                            alpha = alpha_oil;
                            A = A_oil;
                        }
						
                        else {
						
                            abs = stone_fast_abs;
                            scat = stone_fast_scat;
                            ro = ro_stone;
                            alpha = alpha_stone;
                            A = A_stone;
                        }
                    }
                }
			
                r = ((double)rand()/(double)RAND_MAX);
                
                if (r <= scat/(abs + scat)){
				
                    r = ((double)rand()/(double)RAND_MAX);
                    d = cos(r*2*3.141592653);                   //Gives a randomized cos(theta) because z = path length * cos(theta)
                    r = ((double)rand()/(double)RAND_MAX);
                    z = z + (-log(r)/((ro*(abs+scat)*Na/A)*ro*(abs+scat)*Na/A))*d;
                    E = E*((1+alpha)-(1-alpha)*d)/2;
                    interact = interact + 1;
                }
				
                else {
				
                    z = 50;
                    N = N-1;
                    interact = interact + 1;
                    absorb = absorb + 1;
                }
            }
        }while (z <= 40);
    }	

cout<<endl<<"Average number of interactions: "<< interact/No<<endl;
	
cout<<endl<<"Total number past 40cm: "<<N<<endl;
cout<<"From initial number: "<<No<<endl;
cout<<"Number reflected: "<<reflected<<endl;
cout<<"Number absorbed: "<<absorb;
return(0);	
}