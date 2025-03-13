#define _GNU_SOURCE
#include <sched.h>    // for processor affinity
#include <unistd.h>   // unix standard apis
#include <pthread.h>  // pthread api
#include <math.h>

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>

#include "filter.h"
#include "signal.h"
#include "timing.h"

#define MAXWIDTH 40
#define THRESHOLD 2.0
#define ALIENS_LOW  50000.0
#define ALIENS_HIGH 150000.0

int num_threads;
int num_processors;
int fft_size;

int vector_len;       // length of vector we will sum
double* vector;       // the vector we will sum


typedef struct inputs{
  int id;
  double bandwidth;
  int filterOrder;
  double* filterCoeffs;
  signal* sig;
  double* band_power;
  int num_band;
} inputs;

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


void* worker(void* arg) {
  //int blocksize = vector_len / num_threads; // note: floor
  inputs *input = (inputs*)arg;
  int myid = input->id;
  //printf("myid: %d\n", myid);
  

  // put ourselves on the desired processor
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(myid % num_processors, &set);
  if (sched_setaffinity(0, sizeof(set), &set) < 0) { // do it
    perror("Can't setaffinity"); // hopefully doesn't fail
    exit(-1);
  }

  if(myid <= input->num_band - 1){
    generate_band_pass(input->sig->Fs,
                        myid * input->bandwidth + 0.0001, // keep within limits
                        (myid + 1) * input->bandwidth - 0.0001,
                        input->filterOrder,
                        input->filterCoeffs);
    hamming_window(input->filterOrder,input->filterCoeffs);
  
  
    convolve_and_compute_power(input->sig->num_samples,
                              input->sig->data,
                              input->filterOrder,
                              input->filterCoeffs,
                              &(input->band_power[myid]));
  }
  free(input->filterCoeffs); 
  free(input); // prevent memory leak
  // Done.  The master thread will sum up the partial sums
  pthread_exit(NULL);           // finish - no return value
}

unsigned long long int rdtsc(void) {
  unsigned int a;
  unsigned int d;

  // This is inline assembler - rdtsc is the specific instruction
  // we are assembling here
  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));

  return ((unsigned long long)a) | (((unsigned long long)d) << 32);
}

/*
1. remove the dc component from each data point
2. find the average signal power
3. print average signal power to console
4. get the start resources and start the process
5. find start time
6. get the start cycle count
7. make the filter
8. convolution
9. find teh end time and cycle count
10. get the end resources and end the process
11. print the results
*/
int analyze_signal(signal* sig, int filter_order, int num_bands, double* lb, double* ub) {

  double Fc        = (sig->Fs) / 2;
  double bandwidth = Fc / num_bands;

  remove_dc(sig->data,sig->num_samples);

  double signal_power = avg_power(sig->data,sig->num_samples);

  printf("signal average power:     %lf\n", signal_power);

  resources rstart;
  get_resources(&rstart,THIS_PROCESS);
  double time_start = get_seconds();
  unsigned long long tstart = get_cycle_count();

  
  double* band_power = malloc(num_bands * sizeof(double));
  for(int band_index = 0; band_index < num_bands; band_index++){
    band_power[band_index] = -1;
  }
  
  //unsigned long long start = rdtsc();
  double hold_num_bands = (double)num_bands;
  double hold_num_threads = (double)num_threads;
  pthread_t* tid = malloc(num_threads * sizeof(pthread_t));             // array of thread ids
  for (int group = 0; group < ceil(hold_num_bands/hold_num_threads); group++) {
    for (long i = 0; i < num_threads; i++) {
      inputs* thread_inputs = malloc(sizeof(inputs));
      thread_inputs -> id = (group * num_threads)+i;
      thread_inputs->bandwidth = bandwidth;
      thread_inputs->filterOrder = filter_order;
      thread_inputs->filterCoeffs = malloc((filter_order+1)*sizeof(double));
      thread_inputs->sig = sig;
      thread_inputs->band_power = band_power;
      thread_inputs->num_band = num_bands;
      int returncode = pthread_create(&(tid[i]),  // thread id gets put here
                                      NULL, // use default attributes
                                      worker, // thread will begin in this function
                                      (void*)thread_inputs 
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

  }

  unsigned long long tend = get_cycle_count();
  double time_end = get_seconds();

  resources rend;
  get_resources(&rend,THIS_PROCESS);

  resources rdiff;
  get_resources_diff(&rstart, &rend, &rdiff);

  // Pretty print results
  double max_band_power = max_of(band_power,num_bands);
  double avg_band_power = avg_of(band_power,num_bands);
  int wow = 0;
  *lb = -1;
  *ub = -1;

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
        wow = 1;
        if (*lb < 0) {
          *lb = band * bandwidth + 0.0001;
        }
        *ub = (band + 1) * bandwidth - 0.0001;
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
  printf("Analysis took %lf seconds by basic timing\n", time_end - time_start);

  free(tid);
  free(band_power);
  return wow;
}

int main(int argc, char* argv[]) {

  if (argc != 8) {
    usage();
    return -1;
  }

  char sig_type    = toupper(argv[1][0]);
  char* sig_file   = argv[2];
  double Fs        = atof(argv[3]);
  int filter_order = atoi(argv[4]);
  int num_bands    = atoi(argv[5]);
  num_threads = atoi(argv[6]);
  num_processors = atoi(argv[7]);

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

  signal* sig;
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

  sig->Fs = Fs;

  double start = 0;
  double end   = 0;
  if (analyze_signal(sig, filter_order, num_bands, &start, &end)) {
    printf("POSSIBLE ALIENS %lf-%lf HZ (CENTER %lf HZ)\n", start, end, (end + start) / 2.0);
  } else {
    printf("no aliens\n");
  }

  free_signal(sig);

  return 0;
}
