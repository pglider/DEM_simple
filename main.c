#include "decla.h"
#include "functions.c"

int main(){
    
    printf("Parameters : N = %d, R = %f, M = %f, KL = %f, ETA = %f\n", N, R, M, KL, ETA);
    
    //Create output directory
    system("if exist out rd /s /q out");
    system("mkdir out");
    
    //First we initialize the grains
    init();

    //Main time loop
    for(t=1;t<MAX_TIME;t++){
        //Activate the grains
        if((t%ACTIVATION_STEP==0)&&(next_active_grain<N)){
            printf("-------------------------------Adding grain %d\n",next_active_grain-NB);
            disk[next_active_grain].active=1;
            //Give small random horizontal speed
            dx = 1e-7 * ((double) rand()/RAND_MAX-0.5);
            disk[next_active_grain].xold=disk[next_active_grain].x + dx;
            disk[next_active_grain].yold=disk[next_active_grain].y + 2e-7;
            next_active_grain++;
        }

        for(i=0;i<N;i++){
            if(disk[i].active==1){
                //Update the forces
                disk[i].fx=0;// - 1e-5 * disk[i].dx;
                disk[i].fy=-M*G;// - 1e-5 * disk[i].dy;//Gravity
                disk[i].Mz=0;
                disk[i].n_contacts=0;
                if(disk[i].y<-R){
                    disk[i].active=0;
                }
            }
        }

        //Check for collisions
        collision();
        //Integrate the equations of motion
        newton_second_law();

        if(t%EPS_STEP==0){
            printf("Creating eps file %d (%.2f%%)\n",eps_count,100.*t/MAX_TIME);
            r=createps(eps_count);
            if(r==1){
                eps_count++;
            }
            else{
                return 1;
            }
        }

    }
    
}
