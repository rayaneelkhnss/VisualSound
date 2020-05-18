#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>


#define SWAP(a,b) ctmp=(a); (a)=(b); (b)=ctmp
#define ARRAY_LEN(x) ((int) (sizeof (x) / sizeof (x [0])))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define LOG2 10
//#define M_PI 3.14150
#define BUFFER_LEN 1024 // on travaille sur 1024 points
#define HAUTEUR 24 // identique à l'exemple donné
#define LARGEUR	64 // identique à l'exemple donné
// Taille du tableau spectre de 128 car on travaille à 128 000 bits/seconde (fréquence fréquemment utilisée pour le traitement de sons)
#define SPECTRE_LEN 128

double complex TW[BUFFER_LEN];

// Un SNDFILE ressemble beaucoup à un FICHIER dans la bibliothèque Standard C. Les fonctions sf_open_read et sf_open_write renvoient
// un pointeur SNDFILE lorsqu’elles ouvrent avec succès le fichier spécifié.
SNDFILE *infile, *outfile;

// Un pointeur vers une structure SF_INFO est passé à sf_open_read et sf_open_write qui remplissent cette 
// structure d’informations sur le fichier .
SF_INFO sfinfo;

// permet d'inverser les bits de inp (qui contient numbits bits)
int bitrev(int inp, int numbits){
  int rev=0;
  for (int i=0; i<numbits;i++){
      rev=(rev << 1) | (inp & 1);
      inp >>=1;
    }
  return rev;
}

// Cette fonction partitionne le fichier audio en trames et renvoie le nombre de points dans la trame.
sf_count_t sfx_mix_mono_read_double (SNDFILE * file, double * data, sf_count_t datalen){
  SF_INFO info;
  static double multi_data[2048];
  int k, ch, frames_read;
  sf_count_t dataout = 0;
  sf_command (file, SFC_GET_CURRENT_SF_INFO, &info, sizeof(info));
  if (info.channels == 1)
    	return sf_read_double (file, data, datalen);
  while (dataout < datalen){   
   		int this_read ;
		this_read = MIN(ARRAY_LEN(multi_data) / info.channels, datalen);
		frames_read = sf_readf_double (file, multi_data, this_read);
		if (frames_read == 0)
				break;
		for (k = 0 ; k < frames_read ; k++){       
				double mix = 0.0 ;
				for (ch = 0 ; ch < info.channels ; ch++)
						mix += multi_data [k * info.channels + ch] ;
		        data [dataout + k] = mix / info.channels ;
        } ;
    	dataout += this_read ;
    } ;
  	return dataout ;
}

// Fft de manière récursive
void fftrec(double complex *data, double complex *result, unsigned int size, int log2n){
	double complex ypair[size], yimpair[size], Fimpair[size], Fpair[size];
	int n,k,N2;
	if(size>1){
		N2=size/2;
		for(n=0;n<N2;n++){
			ypair[n] = data[n]+data[n+N2];
			yimpair[n] = (data[n]-data[n+N2])*cexp(-2*I*M_PI*n/size);
		}
		fftrec(ypair,Fpair,N2,log2n);
		fftrec(yimpair,Fimpair,N2,log2n);
		for(n=0;n<N2;n++){
			result[2*n]=Fpair[n];
			result[2*n+1]=Fimpair[n];
		}
	}
	else{
		result[0]=data[0];
		return;
	}

}

// permet de précalculer les Twiddle Factors pour simplifier la FFT récursive par la suite
void twiddle(double complex *TW, unsigned int size){
	double complex phi = cexp(-2*I*M_PI/size);
	TW[0]=1;
	for(int i=1; i<size; i++) {
		TW[i]=TW[i-1]*phi;
    }
}

// FFT de façon itérative
void fftiterTW(double complex *data, unsigned int size, int log2n){
	int j, N2, Bpair, Bimpair, Bp=1, N=size;
	double complex impair, pair, ctmp;

	for(int k=0; k<log2n;k++){
		N2=N/2;
		Bpair=0;
		for(int b=0; b<Bp;b++){
			Bimpair=Bpair+N2;
			for(int n=0;n<N2;n++){
				impair = data[Bpair+n] + data[Bimpair+n];
				pair = (data[Bpair+n] - data[Bimpair+n])*TW[n*size/N];
				data[Bpair+n] = pair;
				data[Bimpair+n] = impair;
			}
			Bpair = Bpair+N;
		}
		Bp=Bp*2;
		N=N2;
	}
	for(int i=0;i<size;i++){
		j=bitrev(i,log2n);
		if(j>i){
			SWAP(data[j],data[i]);
		}
	}
	for(int i=size-1;i>0;i--){
		data[i]=data[i-1];
	}
	data[0]=ctmp;
	return;
}

// DFT classique d'un signal 
void DFT(double complex *data, double complex *dft, unsigned int size){
	for(int i=0; i<BUFFER_LEN;i++){
		for(int k=0; k<BUFFER_LEN;k++){
			dft[i]+=data[k]*cexp(-2*I*M_PI*k*i/size);
		}
	}
}

// convertit le nombre en complexe
double complex doubleToComplexe(double nombre){
	return nombre+I*0.0;
}

// renvoie le temps d'éxécution d'un sous programme 
double timing(double start, double stop) {
    return (double)(stop-start)/(double)(CLOCKS_PER_SEC);
}

// renvoie le max d'un tableau à k éléments
double max(double t[], double k){
	double max=t[0]; // on initialise le max au premier élément
	for(int i=1;i<k;i++){
		// si l'élément parcouru est isupérieur au max alors il devient le nouveau max
		if(t[i]>max){
			max=t[i];
		}
	}
	return max;
}

// renvoie le module d'un complexe
double module(double complex c){
    return sqrt(creal(c)*creal(c) + cimag(c)*cimag(c));
}

// permet de définir la durée entre chaque affichage du spectre
void stop(int t){
	clock_t begin=clock();
	while((clock()-begin) <= (t*CLOCKS_PER_SEC/1000)); //temps écoulé en réel différent du nombre de ticks consommé (programme peut passer une partie du temps à ne rien faire)
}

// on créé le spectre
void initSpectre(double complex dataComplex[],double spectre[],int size){
	int k=0;
	int abs=0;  	
	double amplitude=0.0; 
	// on parcourt le tableau dataComplex
	for(int i=0; i<size;i++){
		amplitude+=module(dataComplex[i]); // l'amplitude pour un point donné est la somme des modules des valeurs complexes du tableau
		k++;
		if(k==size/SPECTRE_LEN){ // lorsque l'on atteint le dernier point de la trame on a notre amplitude finale pour la trame
			spectre[abs]=20*log10(amplitude/k); // on affecte l'amplitude au spectre pour le point d'abscisse abs
			k=amplitude=0; // on réinitialise k et amplitude à 0
			abs++; // on travaillera ensuite sur la prochaine abscisse
		}
	}
}

// on affiche le spectre à l'écran
void printSpectre(double *spectre){
	double maximum=max(spectre,LARGEUR);		
	double ecran[LARGEUR];
	for(int i=0;i<LARGEUR;i++){
		ecran[i]=spectre[i];
	}
	for(int i=0;i<HAUTEUR;i++){
		for(int k=0;k<LARGEUR;k++){
			if(ecran[k]>=maximum-i){ // si la valeur du spectre en un point est supérieure ou égale au max des valeurs prises par le spectre
				printf("*");        // alors on la représente (*), l'étoile la plus haute ne pourra être atteinte que là où la valeur du spectre
			}                        // est maximale. Plus on sera vers le bas de l'écran plus i sera grand et ainsi il sera plus fréquent de
			else  printf(" ");       // majorer (maximum-i) , ce qui est logique.
		}
		printf("\n");
	}	
}


int affichage_fft(){
    double data[BUFFER_LEN];
	// on traitera les données de type complexe (conversion préalable)
	double complex dataComplex[BUFFER_LEN];
	double spectre[SPECTRE_LEN];
	sf_count_t readcount;
	twiddle(TW, BUFFER_LEN); //on exécute twiddle pour alléger les futurs calculs de la fft
	infile = sf_open("ordinary_speech.wav",SFM_READ, &sfinfo);
	// son impossible à ouvrir
	if(infile==NULL){
		fprintf(stderr, "Fichier d'entrée impossible à ouvrir\n");
		return EXIT_FAILURE;
	}
	double channels = sfinfo.channels;
    double samplerate = sfinfo.samplerate;
    // tant qu'il y a des trames on les lit puis on fait leur FFT pour avoir le spectre que l'on affichera à l'écran
    // un nombre de points dans la trame nul indiquera que l'on a finit de lire le fichier
	while((readcount=sfx_mix_mono_read_double(infile,data,BUFFER_LEN))>0){
		//conversion du tableau de doubles en un tableau de complexes
		for(int i=0; i<BUFFER_LEN; i++){
			dataComplex[i] = doubleToComplexe(data[i]);
		}
		fftiterTW(dataComplex,BUFFER_LEN,LOG2); // on utilise la FFT itérative sur chaque trame
		initSpectre(dataComplex,spectre,BUFFER_LEN); //le spectre est créé
		system("clear");
		printSpectre(spectre); //on affiche sur la sortie standard du terminal le spectre
		stop(120); //durée entre chaque affichage du spectre
	}
	sf_close(infile); //on ferme le fichier d'entrée
	return EXIT_SUCCESS;
}

void temps_de_calcul(){
	// on ouvre le fichier d'entrée
	infile = sf_open("Do-sol.wav", SFM_READ, &sfinfo);
    if (infile == NULL) {
        printf("Fichier d'entrée impossible à ouvrir \n");
        sf_perror(NULL);
    }
    // on ouvre le fichier de sortie
    outfile = sf_open("pitch_in_air.wav",SFM_WRITE, &sfinfo);
    if (outfile == NULL) {
        printf("Fichier de sortie impossible à ouvrir \n");
        sf_perror(NULL);
    }
    double data[BUFFER_LEN];
    // on traitera les données de type complexe (conversion préalable)
	double complex dataComplex[BUFFER_LEN];
	double channels = sfinfo.channels;
    double samplerate = sfinfo.samplerate;
    sf_count_t readcount;
	double tempsDFT, tempsFFTiterTW, tempsFFTfrec;
	double complex dft[BUFFER_LEN];
	double complex spectre[BUFFER_LEN];
   
   // on va calculer les temps d'éxécution sur la première trame car on devrait obtenir les mêmes temps d'exécution peu importe la trame choisie
    readcount = sfx_mix_mono_read_double(infile, data, BUFFER_LEN); 
    for(int i=0; i<BUFFER_LEN;i++){ 
        dataComplex[i] = doubleToComplexe(data[i]); 
    };
 
	//Temps d'éxécution pour la FFT itérative
	struct timeval temps_avant_1, temps_apres_1; //On définit deux varibales pour le temps avant et après qui sont des pointeurs donc elles contiennent 2 champs
    gettimeofday (&temps_avant_1, NULL); //on prend le temps avant le calcul de la fft
    twiddle(TW,BUFFER_LEN);//on exécute twiddle pour alléger les futurs calculs de la fft
	fftiterTW(dataComplex,BUFFER_LEN,LOG2);//on exécute la fft avec twiddle
    gettimeofday (&temps_apres_1, NULL);//on prend la temps tout de suite après la fin de la ftttw
    tempsFFTiterTW=timing(temps_avant_1.tv_usec, temps_apres_1.tv_usec);//on fait le temps stop-le temps de start et on a le temps d'exècution !
	//La date fournie par la fonction "gettimeofday" est stockée dans une structure de type "struct timeval" qui contient les champs tv_sec et tv_usec correspondants au temps courant récupéré.
	//Donc nous prenons .tv_usec puisque qu'il s'agit des microsecondes
    	
	//Temps d'éxécution pour la FFt récursive
	struct timeval temps_avant_2, temps_apres_2;
    gettimeofday (&temps_avant_2, NULL);  
    fftrec(dataComplex,spectre,BUFFER_LEN,LOG2);//on exécute la fft récursive
    gettimeofday (&temps_apres_2, NULL);   	
    tempsFFTfrec=timing(temps_avant_2.tv_usec, temps_apres_2.tv_usec);//De même      

	//Temps d'exécution pour une DFT classique
	struct timeval temps_avant_3, temps_apres_3;
	gettimeofday(&temps_avant_3, NULL);
	DFT(dataComplex, dft, BUFFER_LEN);//on exécute la DFT classique 
	gettimeofday(&temps_apres_3, NULL);
	tempsDFT = timing(temps_avant_3.tv_usec, temps_apres_3.tv_usec);//De même
        
    
	// fermeture des fichiers d'entrée et de sortie
    sf_close(infile);
    sf_close(outfile);

    printf("Durée d'une FFT itérative : %lf ms\n",tempsFFTiterTW*1000);
    printf("Durée d'une FFT récursive : %lf ms\n",tempsFFTfrec*1000);
    printf("Durée d'une DFT : %lf ms\n", tempsDFT*1000);

}



int main(int argc, char*argv[]){
	affichage_fft();
	temps_de_calcul();
	
}


