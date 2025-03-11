#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#include <sched.h>    // for processor affinity
#include <unistd.h>   // unix standard apis
#include <pthread.h>  // pthread api
#include <math.h>

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0


double* filterCoeffs;
int filter_order;
int num_bands;

signal* sig;                // signal 
int vector_len;       // length of the signal vector

int component_width;
int component_vector_length;
double** component_vector;       // the vector that will contain the components of the filtered signal
double* output_vector;
double** combine_vector;

int num_threads;            // number of threads we will use
int num_processors;         // number of processors we will use
pthread_t* tid;             // array of thread ids

void usage() {
  printf("usage: band_scan text|bin|mmap signal_file Fs filter_order num_bands\n");
}


double avg_power(double* data, int num) {

  double ss = 0;
  for (int i = 0; i < num; i++) {
    ss += data[i] * data[i];
  }

  return ss / num;
}

double max_of(double* data, int num) {

  double m = data[0];
  for (int i = 1; i < num; i++) {
    if (data[i] > m) {
      m = data[i];
    }
  }
  return m;
}

double avg_of(double* data, int num) {

  double s = 0;
  for (int i = 0; i < num; i++) {
    s += data[i];
  }
  return s / num;
}

void remove_dc(double* data, int num) {

  double dc = avg_of(data,num);

  printf("Removing DC component of %lf\n",dc);

  for (int i = 0; i < num; i++) {
    data[i] -= dc;
  }
}

// Function run by each thread
void* worker(void* arg) {
    long myid     = (long)arg;
    // put ourselves on the desired processor
    cpu_set_t set;
    CPU_ZERO(&set);
    CPU_SET(myid % num_processors, &set);
    if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
      perror("Can't setaffinity"); // hopefully doesn't fail
      exit(-1);
    }
    // This figures out the chunk of the vector I should
  // work on based on my id
  int mystart = myid * component_width;
  int myend   = myid*(component_width + 1);

  double* output = {component_width * 2};
  printf("worker pre\n");
  fft_convolute(component_width, 
               ((sig->data) + mystart), 
               filter_order, 
               filterCoeffs,
               output);
  printf("worker post\n");

  component_vector[myid] = output;

  // Done.  The master thread will sum up the partial sums
  pthread_exit(NULL);           // finish - no return value

}

int main(int argc, char* argv[]) {
/////////////////////////////////////// Read inputs//////////////////////////////////////
    if (argc != 8) {
        usage();
        return -1;
    }
    
    char sig_type    = toupper(argv[1][0]);
    char* sig_file   = argv[2];
    double Fs        = atof(argv[3]);
    filter_order = atoi(argv[4]);
    num_bands    = atoi(argv[5]);
    num_threads  = atoi(argv[6]);            // number of threads we will use
    num_processors = atoi(argv[7]);           // number of processors we will use
    
    assert(Fs > 0.0);
    assert(filter_order > 0 && !(filter_order & 0x1));
    assert(num_bands > 0);
    
    printf("type:     %s\n\
    file:     %s\n\
    Fs:       %lf Hz\n\
    order:    %d\n\
    bands:    %d\n",
             sig_type == 'T' ? "Text" : (sig_type == 'B' ? "Binary" : (sig_type == 'M' ? "Mapped Binary" : "UNKNOWN TYPE")),
             sig_file,
             Fs,
             filter_order,
             num_bands);
    printf("Load or map file\n");
///////////////////////////////////// End Read Inputs //////////////////////////////////////////////////////////////////
///////////////////////////////////// Load the correct signal //////////////////////////////////////////////////////////

    switch (sig_type) {
      case 'T':
        sig = load_text_format_signal(sig_file);
        break;
  
      case 'B':
        sig = load_binary_format_signal(sig_file);
        break;
  
      case 'M':
        sig = map_binary_format_signal(sig_file);
        break;
  
      default:
        printf("Unknown signal type\n");
        return -1;
    }
  
    if (!sig) {
      printf("Unable to load or map file\n");
      return -1;
    }
//////////////////////////////////////////// End Load Correct Signal ///////////////////////////////////////////////////////
double signal_power = avg_power(sig->data, sig->num_samples);
printf("signal average power:     %lf\n", signal_power);
component_width = (sig->num_samples) / num_threads;
int fft_size = 1;
while(fft_size < component_width){ // find correct component width
  fft_size *= 2;
}
component_width = fft_size;
component_vector_length = ceil((sig->num_samples) / component_width);
int last_ele_input_length = sig->num_samples - (component_width * (component_vector_length - 1));

//////////////////////////////////////////// Memory Allocation for threads /////////////////////////////////////////////////
    vector_len  = sig->num_samples;    // length of signal vector 
  
    component_vector      = (double*)malloc(sizeof(double) * vector_len);
    tid         = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);
  
    if (!component_vector || !tid) {
      fprintf(stderr, "cannot allocate memory\n");
      exit(-1);
    }
/////////////////////////////////////////////End Memory Allocation for threads///////////////////////////////////////////////
///////////////////////////////////////////// Thread Creation and Run //////////////////////////////////////////////////////
  
  sig->Fs = Fs;
  double Fc        = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;
  filterCoeffs[filter_order + 1];
  double band_power[num_bands];
  combine_vector[num_bands];

  remove_dc(sig->data,sig->num_samples);
  
  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  

  for(int band = 0; band < num_bands; band++) {
    printf("gen bandpass\n");
    generate_band_pass(sig->Fs,
      band * bandwidth + 0.0001, // keep within limits
      (band + 1) * bandwidth - 0.0001,
      filter_order,
      filterCoeffs);
    printf("pre hamming\n");
    hamming_window(filter_order,filterCoeffs);
    printf("post hamming\n");
    // launch threads
    // thread i will compute partial_sum[i], which will sum
    // vector[i*ceiling(vector_size/num_threads)
    //          through
    //        (i+1)*floor(vector_size/num_threads) ]
    // the last thread will also handle the additional elements
    for (long i = 0; i < num_threads; i++) {
      int returncode = pthread_create(&(tid[i]),  // thread id gets put here
                                      NULL, // use default attributes
                                      worker, // thread will begin in this function
                                      (void*)i // we'll give it i as the argument
                                      );
      if (returncode != 0) {
        perror("Failed to start thread");
        exit(-1);
      }
    }
    // now we will join all the threads
    for (int i = 0; i < num_threads; i++) {
      int returncode = pthread_join(tid[i], NULL);
      if (returncode != 0) {
        perror("join failed");
        exit(-1);
      }
    }

    for(int ii = 0; ii < component_vector_length-1; ii++){ //all but last element
      for(int jj = 0; jj < sizeof(component_vector[ii]); jj++){
        output_vector[(ii*component_width)+jj] += component_vector[ii][jj];
      }
    }

    int last_ele_length = sizeof(component_vector[component_vector_length-1]);
    int counter = 0;
    for(int ii = sizeof(output_vector) - (last_ele_input_length + filter_order); ii < sizeof(output_vector); ii++){
      output_vector[ii] += component_vector[component_width - 1][counter];
      counter++;
    }

    combine_vector[band] = output_vector;
    band_power[band] = avg_power(output_vector, sizeof(output_vector));
  }
  unsigned long long tend = get_cycle_count();
  double end = get_seconds();
  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);
  ////////////////////////////////////////////// End Thread Creation and Run //////////////////////////////////////
  ////////////////////////////////////////////// Combine Thread Results ///////////////////////////////////////////
  
  /////////////////////////////////////////////End Combine Thread Results ////////////////////////////////////////
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int found = 0;
  double lb = -1;
  double ub   = -1;
  
  for (int band = 0; band < num_bands; band++) {
    double band_low  = band * bandwidth + 0.0001;
    double band_high = (band + 1) * bandwidth - 0.0001;

    printf("%5d %20lf to %20lf Hz: %20lf ",
           band, band_low, band_high, band_power[band]);

    for (int i = 0; i < MAXWIDTH * (band_power[band] / max_band_power); i++) {
      printf("*");
    }

    if ((band_low >= ALIENS_LOW && band_low <= ALIENS_HIGH) ||
        (band_high >= ALIENS_LOW && band_high <= ALIENS_HIGH)) {

      // band of interest
      if (band_power[band] > THRESHOLD * avg_band_power) {
        printf("(WOW)");
        found = 1;
        if (lb < 0) {
          lb = band * bandwidth + 0.0001;
        }
        ub = (band + 1) * bandwidth - 0.0001;
      } else {
        printf("(meh)");
      }
    } else {
      printf("(meh)");
    }

    printf("\n");
  }

  printf("Resource usages:\n\
  User time        %lf seconds\n\
  System time      %lf seconds\n\
  Page faults      %ld\n\
  Page swaps       %ld\n\
  Blocks of I/O    %ld\n\
  Signals caught   %ld\n\
  Context switches %ld\n",
         rdiff.usertime,
         rdiff.systime,
         rdiff.pagefaults,
         rdiff.pageswaps,
         rdiff.ioblocks,
         rdiff.sigs,
         rdiff.contextswitches);

  printf("Analysis took %llu cycles (%lf seconds) by cycle count, timing overhead=%llu cycles\n"
         "Note that cycle count only makes sense if the thread stayed on one core\n",
         tend - tstart, cycles_to_seconds(tend - tstart), timing_overhead());
  printf("Analysis took %lf seconds by basic timing\n", end - start);

  if(found == 1){
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", lb, ub, (ub + lb) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);
  
    return 0;
}
  
  