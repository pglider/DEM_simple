#include "decla.h"

void init(){
    int i;
    double r_poly;
    //Bottom grains
    for(i=0;i<NB;i++){
        disk[i].x = i * 2 * R * 1.1;
        disk[i].y = 0.;
        disk[i].Oz = 0.;
        disk[i].r = R;
        disk[i].mass = M;
        disk[i].inertia = 2./5 * M * R * R;
        disk[i].fixed=1;
        disk[i].active=1;
        disk[i].yold=disk[i].y;
        disk[i].xold=disk[i].x;
        disk[i].Ozold=disk[i].Oz;
        disk[i].dx=0;
        disk[i].dy=0;
        disk[i].dOz=0;
    }
    //Rest of the grains
    for(i=NB;i<N;i++){
        disk[i].x = (NB-1) * R * 1.1;//Middle of the bottom grains
        disk[i].y = Y0 * R;//Initial height
        disk[i].Oz = 0.;
        r_poly = 1. + 0.2 * ((double) rand()/RAND_MAX-0.5);//20% polydispersity
        disk[i].r = R * r_poly;
        disk[i].mass = M * r_poly * r_poly * r_poly;
        disk[i].inertia = 2./5 * M * R * R * r_poly * r_poly * r_poly * r_poly * r_poly;
        disk[i].fixed=0;
        disk[i].active=0;
        disk[i].yold=disk[i].y;
        disk[i].xold=disk[i].x;
        disk[i].Ozold=disk[i].Oz;
        disk[i].dx=0;
        disk[i].dy=0;
        disk[i].dOz=0;
    }
    next_active_grain=NB;

}

void collision(){
    int i, j, current_contact, contact_index;
    double dx, dy, distance, delta, fn, nx, ny, vn, vt, vtx, vty, ftx ,fty;
    double utx, uty, ut, vtx_n, vty_n;

    for (i = 0; i < N; i++) {
        if (disk[i].active == 1) {
            for (j = i + 1; j < N; j++) {
                if (disk[j].active == 1) {

                    //First determine if the grains are in contact
                    dx = disk[j].x - disk[i].x;
                    dy = disk[j].y - disk[i].y;
                    distance = sqrt(dx * dx + dy * dy);
                    delta = disk[i].r + disk[j].r - distance;

                    if (delta > 0) {
                        // Collision detected
                        // Looking if the contact is already existing
                        for (contact_index=0;contact_index<MAXCONTACTS;contact_index++) {
                            if (disk[i].contactsold[contact_index]==j) {//Contact already existing : load previous values of ut
                                disk[i].utx[disk[i].n_contacts]=disk[i].utxold[contact_index];
                                disk[i].uty[disk[i].n_contacts]=disk[i].utyold[contact_index];
                                break;
                            }
                        }
                        current_contact = disk[i].n_contacts;
                        disk[i].contacts[current_contact]=j;
                        disk[i].n_contacts++;

                        // Definition of the unit vector normal to the contact
                        nx=dx/distance;
                        ny=dy/distance;

                        // Compute the normal speed Vn=\vec{v}.\vec{n} for the viscoelasic model
                        // = projection of normal speed in the direction of the normal vector
                        vn=(disk[j].dx-disk[i].dx)*nx + (disk[j].dy-disk[i].dy)*ny;
                        
                        //Normal force : viscoelastic model
                        fn = - (KL* delta - ETA*vn);

                        // Compute the tangential speed vt using Varignon's formula
                        vtx = disk[j].dx - disk[i].dx - vn*nx + ny*(disk[j].dOz*disk[j].r + disk[i].dOz*disk[i].r);
                        vty = disk[j].dy - disk[i].dy - vn*ny - nx*(disk[j].dOz*disk[j].r + disk[i].dOz*disk[i].r);
                        vt= sqrt(vtx*vtx+vty*vty);

                        // Cundall's model
                        // Compute the integral displacement ut, here expressed in N
                        utx=disk[i].utx[current_contact]+KT*vtx*DT;
                        uty=disk[i].uty[current_contact]+KT*vty*DT;
                        ut=sqrt(utx*utx+uty*uty);

                        //Test if we are in the sliding regime or not
                        if(fabs(ut/(MU*fn))>1. && (utx*vtx)>0 && (uty*vty)>0) {//Sliding regime
                            //First we determine the direction of the normalized sliding speed
                            if(vt==0.) {
                                vtx_n=0.;
                                vty_n=0.;
                            }
                            else {
                                vtx_n=vtx/vt;
                                vty_n=vty/vt;
                            }
                            // Then we determine the friction force
                            ftx=- MU*fn*vtx_n;
                            fty=- MU*fn*vty_n;

                            disk[i].utx[current_contact]=ftx;
                            disk[i].uty[current_contact]=fty;
                        }
                        else{//Non sliding regime
                            disk[i].utx[current_contact]+=KT*vtx*DT;
                            disk[i].uty[current_contact]+=KT*vty*DT;
                            ftx=disk[i].utx[current_contact];// + 0.001*vtx;
                            fty=disk[i].uty[current_contact];// + 0.001*vty;
                        }

                        // Update forces and moments

                        disk[i].fx+=fn*nx+ftx;
                        disk[i].fy+=fn*ny+fty;
                        disk[i].Mz+= disk[i].r * (nx*fty-ny*ftx);

                        disk[j].fx+=-fn*nx-ftx;
                        disk[j].fy+=-fn*ny-fty;
                        disk[j].Mz+= disk[j].r * (nx*fty-ny*ftx);

                    }
                }
            }
        }
    }
}

void newton_second_law(){
    int i;
    double new_x, new_y, new_Oz;
    for(i=0;i<N;i++){
        if(disk[i].active==1){
            if(disk[i].fixed==1){
                disk[i].xold = disk[i].x;
                disk[i].yold = disk[i].y;
                disk[i].dx=0;
                disk[i].dy=0;
            }
            else{
                // Update positions using Verlet integration
                new_x = 2 * disk[i].x - disk[i].xold + (disk[i].fx / disk[i].mass) * DT * DT;
                new_y = 2 * disk[i].y - disk[i].yold + (disk[i].fy / disk[i].mass) * DT * DT;
                new_Oz = 2 * disk[i].Oz - disk[i].Ozold + (disk[i].Mz / disk[i].inertia) * DT * DT;

                // Update previous positions
                disk[i].xold = disk[i].x;
                disk[i].yold = disk[i].y;
                disk[i].Ozold = disk[i].Oz;

                //Update speeds
                disk[i].dx=(new_x-disk[i].xold)/(2*DT);
                disk[i].dy=(new_y-disk[i].yold)/(2*DT);
                disk[i].dOz=(new_Oz-disk[i].Ozold)/(2*DT);

                // Update current positions
                disk[i].x = new_x;
                disk[i].y = new_y;
                disk[i].Oz = new_Oz;
            }
        }
        else{
            disk[i].xold = disk[i].x;
            disk[i].yold = disk[i].y;
            disk[i].Ozold = disk[i].Oz;
            disk[i].dx=0;
            disk[i].dy=0;
            disk[i].dOz=0;

        }
    }
}

int createps(int ieps){
    char name[255];// Name of the eps file
    double xref,yref;// Reference position for the grains
    int i;// Counter for the grains
    sprintf(name,"out\\p%.5d.eps",ieps);
    eps_file=fopen(name,"w");
    if(eps_file==NULL){
        printf("Error opening file %s\n",name);
        return -1;
    }

    // Header
    fprintf(eps_file,"%%!PS-Adobe-3.0 EPSF-3.0\n");
    fprintf(eps_file,"%%");
    fprintf(eps_file,"%%BoundingBox: 0 0 %d %d\n", (int)(2*(NB-1)*1.1*R*SCALE), (int)((Y0+2)*R*SCALE));
    fprintf(eps_file,"%%");
    fprintf(eps_file,"%%EndComments \n");

    // Shortcuts
    fprintf(eps_file,"/M{moveto}def\n");				// Move to
    fprintf(eps_file,"/L{lineto stroke}def\n");			// Line to
    fprintf(eps_file,"/D{1 1 sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore fill stroke}def\n");//Draw disk
    fprintf(eps_file,"/G{0 0 0.5 sethsbcolor newpath 0 360 arc closepath gsave 0 setgray grestore fill stroke}def\n");//Draw gray disk
    fprintf(eps_file,"/C{0 0 0 sethsbcolor newpath 0 360 arc closepath stroke}def\n");//Draw circle
    
    fprintf(eps_file,"1 setlinewidth\n");
    
    for(i=0;i<=N;i++){		
        xref=disk[i].x;
        yref=disk[i].y;
        if(disk[i].fixed==1){
            fprintf(eps_file,"\n %f %f %f %f D", xref*SCALE, yref*SCALE, disk[i].r*SCALE, 0.4);            
            fprintf(eps_file,"\n %f %f %f C", xref*SCALE, yref*SCALE, disk[i].r*SCALE);
        }
        else{
            if(disk[i].active==1){
                fprintf(eps_file,"\n %f %f %f %f D", xref*SCALE, yref*SCALE, disk[i].r*SCALE, 0.8);
                fprintf(eps_file,"\n %f %f %f C", xref*SCALE, yref*SCALE, disk[i].r*SCALE);
                // Draw line showing orientation
                fprintf(eps_file,"\n %f %f M %f %f L", xref*SCALE, yref*SCALE,(xref + disk[i].r * cos(disk[i].Oz))*SCALE, (yref + disk[i].r * sin(disk[i].Oz))*SCALE);
            }
            else{
                fprintf(eps_file,"\n %f %f %f G", xref*SCALE, yref*SCALE, disk[i].r*SCALE);
            }
        }
    }
    fclose(eps_file);
    return 1;
}