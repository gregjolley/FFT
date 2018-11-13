#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define size  16384
#define log_size 14
#define PI  3.1415926535897932384626433832795

// FFT code development with intended porting to PIC/dsPIC24 devices, or something more capable


void separate(double signal[]){

    double temp[size];
    int i, j, k, l, power;

    for(i=1; i<log_size; i++){
        power = 1;
        for(l=0; l<i; l++){
            power = power*2;    // 2, 4, 8,...
        }
        for(j=0; j<(power/2); j++){ // 1, 2, 4,...cycles
            for(k=0; k<(size/power); k++){ // 1/2, 1/4, 1/8,...
                temp[k + j*2*size/power] = signal[2*k + j*2*size/power];
                temp[k + (j*2+1)*size/power] = signal[2*k+1 + j*2*size/power];
            }
        }
        for(l=0; l<size; l++){
            signal[l] = temp[l];
        }
    }
}


void FFT_1(double signal[], double FFT[]){

    int i, k, l;
    double temp_Re[size], temp_Im[size];
    double signal_Re[size], signal_Im[size];


    separate(signal);


    // ********** Calculate the 1st stage N = 2 **********

    for(i=0; i<size/2; i++){
        temp_Im[i*2] = temp_Im[i*2+1] = 0.0;
        temp_Re[i*2] = signal[i*2] + signal[i*2+1];
        temp_Re[i*2+1] = signal[i*2] - signal[i*2+1];
    }

    for(i=0; i<size; i++){
        signal_Re[i] = temp_Re[i];
        signal_Im[i] = temp_Im[i];
    }


    for(l=4; l<=size; l=l*2){            // 4, 8, 16,...,size
        for(i=0; i<(size/l); i++){
            for(k=0; k<(l/2); k++){         // 0, 1,...,l/2-1

                temp_Re[i*l + k] = signal_Re[i*l + k] + signal_Re[i*l + k+l/2]*cos(2.0*PI*((double)k)/((double)l));     // E + O*w

                temp_Re[i*l + k] -= signal_Im[i*l + k+l/2]*sin(-2.0*PI*((double)k)/((double)l));        // O*w

                temp_Im[i*l + k] = signal_Im[i*l + k] + signal_Im[i*l + k+l/2]*cos(2.0*PI*((double)k)/((double)l));     // E + O*w

                temp_Im[i*l + k] += signal_Re[i*l + k+l/2]*sin(-2.0*PI*((double)k)/((double)l));    // O*w

                temp_Re[i*l + k+l/2] = signal_Re[i*l + k] - signal_Re[i*l + k+l/2]*cos(2.0*PI*((double)k)/((double)l));     // E - O*w

                temp_Re[i*l + k+l/2] += signal_Im[i*l + k+l/2]*sin(-2.0*PI*((double)k)/((double)l));        // -O*w

                temp_Im[i*l + k+l/2] = signal_Im[i*l + k] - signal_Im[i*l + k+l/2]*cos(2.0*PI*((double)k)/((double)l));     // E - O*w

                temp_Im[i*l + k+l/2] -= signal_Re[i*l + k+l/2]*sin(-2.0*PI*((double)k)/((double)l));        // -O*w

            }
        }
        for(i=0; i<size; i++){
            signal_Re[i] = temp_Re[i];
            signal_Im[i] = temp_Im[i];
        }
    }

     for(i=0; i<size; i++){
        FFT[i] = sqrt(signal_Re[i]*signal_Re[i] + signal_Im[i]*signal_Im[i]);
    }

}




void FFT_2(double signal[], double FFT[]){

    int i, k, l;
    double temp_Re1, temp_Im1, temp_Re2, temp_Im2;
    double signal_Re[size], signal_Im[size];
    double temp_cos, temp_sin;


    separate(signal);


    // ********** Calculate the 1st stage N = 2 **********

    for(i=0; i<size/2; i++){

        temp_Re1 = signal[i*2] + signal[i*2+1];
        temp_Re2 = signal[i*2] - signal[i*2+1];

        signal[i*2] = temp_Re1;
        signal[i*2+1] = temp_Re2;
    }

    for(i=0; i<size; i++){
        signal_Re[i] = signal[i];
        signal_Im[i] = 0.0;
    }


    for(l=4; l<=size; l=l*2){            // 4, 8, 16,...,size
        for(i=0; i<(size/l); i++){
            for(k=0; k<(l/2); k++){         // 0, 1,...,l/2-1

                temp_cos = cos(2.0*PI*((double)k)/((double)l));
                temp_sin = sin(-2.0*PI*((double)k)/((double)l));

                temp_Re1 = signal_Re[i*l + k] + signal_Re[i*l + k+l/2]*temp_cos;     // E + O*w

                temp_Re1 -= signal_Im[i*l + k+l/2]*temp_sin;        // O*w

                temp_Im1 = signal_Im[i*l + k] + signal_Im[i*l + k+l/2]*temp_cos;     // E + O*w

                temp_Im1 += signal_Re[i*l + k+l/2]*temp_sin;    // O*w

                temp_Re2 = signal_Re[i*l + k] - signal_Re[i*l + k+l/2]*temp_cos;     // E - O*w

                temp_Re2 += signal_Im[i*l + k+l/2]*temp_sin;        // -O*w

                temp_Im2 = signal_Im[i*l + k] - signal_Im[i*l + k+l/2]*temp_cos;     // E - O*w

                temp_Im2 -= signal_Re[i*l + k+l/2]*temp_sin;        // -O*w

                signal_Re[i*l + k] = temp_Re1;
                signal_Re[i*l + k+l/2]= temp_Re2;
                signal_Im[i*l + k] = temp_Im1;
                signal_Im[i*l + k+l/2]= temp_Im2;


            }
        }
    }

     for(i=0; i<size; i++){
        FFT[i] = sqrt(signal_Re[i]*signal_Re[i] + signal_Im[i]*signal_Im[i]);
    }

}



void DFT(double signal[], double FFT[]){

    int i, k;
    double temp_Re[size], temp_Im[size];

    for(k=0; k<size; k++){
        temp_Re[k] = temp_Im[k] = 0.0;
        for(i=0; i<size; i++){
            temp_Re[k] += signal[i]*cos(2.0*PI*(double)(i*k)/(double)size);
            temp_Im[k] -= signal[i]*sin(2.0*PI*(double)(i*k)/(double)size);
        }
    }
    for(i=0; i<size; i++){
        FFT[i] = sqrt(temp_Re[i]*temp_Re[i] + temp_Im[i]*temp_Im[i]);
    }

}


int main(){

    int i, j, temp, k, l;
    double signal[size];    // load with some input data
    double FFT[size], Re[size], Im[size];
    double temp_Re[size], temp_Im[size];

    double signal_Re[size], signal_Im[size];

    // Load array with windowed sine signal
    for(i=0; i<size; i++){
        signal[i] = 0.5*(1-cos((2.0*PI*(double)i/(double)(size-1))))*sin(2.0*PI*15.0*(double)i/(double)size);
        //signal[i] = sin(2.0*PI*51.5*(double)i/(double)size);
    }

    //for(i=0; i<size; i++){
    //    signal[i] = i;
    //}


    clock_t begin = clock();
    DFT(signal, FFT);
    clock_t end = clock();

    printf("%f\n", (double)(end-begin));
    for(l=0; l<30; l++){
        printf("%d  %f\n", l, FFT[l]);
    }

    begin = clock();
    FFT_2(signal, FFT);
    end = clock();

    printf("%f\n", (double)(end-begin));

    for(l=0; l<30; l++){
        printf("%d  %f\n", l, FFT[l]);
    }
}
