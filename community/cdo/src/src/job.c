#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600 /* gethostname */
#endif

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>

#include <sys/types.h>  /* fstat */
#include <sys/stat.h>
#include <unistd.h>

#include "cdo.h"

#if defined(HAVE_LIBDRMAA)
#  include "drmaa.h"
#endif

#define  GRID_TMPDIR  "/opt/griddata/tmp"

int ftpget(int flag, const char *url, const char *path, const char *target, const char *source);

#if defined(HAVE_LIBDRMAA)
static drmaa_job_template_t *create_job_template(const char *expname, const char *jobfilename, const char *jobname, const char *tmppath)
{
  drmaa_job_template_t *job = NULL;

  char error[DRMAA_ERROR_STRING_BUFFER];
  char name[DRMAA_ATTR_BUFFER], value[DRMAA_ATTR_BUFFER];
  char attr[1024];

  long size;
  char *dir, *ptr;

  drmaa_attr_names_t *job_attributes;

  int drmaa_errno;

  char host[1024];

  char *output_path;

  int len, len1, len2;

  /* determine hostname */

#if defined(HAVE_GETHOSTNAME)
  gethostname(host, sizeof(host));  
#else
  fprintf(stderr, "Function gethostname not available!\n");
  exit(EXIT_FAILURE);
#endif

  /* determine current path */

  size = pathconf(".", _PC_PATH_MAX);
  if ( (dir = (char*) malloc((size_t)size)) != NULL )
    {
      ptr = getcwd(dir, (size_t)size);
    }

  /* generate DRMAA conform output path */
  
  len1 = strlen(host);
  /*len2 = strlen(dir);*/
  len2 = strlen(GRID_TMPDIR);
  len = len1+len2+2;

  output_path = (char*) malloc(len*sizeof(char));
  /*
  strcpy(output_path, host);
  strcat(output_path, ":");
  */
  strcpy(output_path, ":");
  strcat(output_path, GRID_TMPDIR);

  /* need to allow chdir on execution host, not thread save! */

  /* setenv("SGE_DRMAA_ALLOW_CWD", "yes", 1); */

  /* allocate job template */

  if (drmaa_allocate_job_template(&job, NULL, 0) != DRMAA_ERRNO_SUCCESS)
    return NULL;

  /* the job's name */
  drmaa_set_attribute(job, DRMAA_JOB_NAME, jobname, NULL, 0);

  /* the job to be run */
  drmaa_set_attribute(job, DRMAA_REMOTE_COMMAND, jobfilename, NULL, 0);

  /* submit state */
  drmaa_set_attribute(job, DRMAA_JS_STATE, "drmaa_active", NULL, 0);

  /* working directory on execution host */
  /* drmaa_set_attribute(job, DRMAA_WD, GRID_TMPDIR, NULL, 0); */
  drmaa_set_attribute(job, DRMAA_WD, tmppath, NULL, 0);

  /* path for output */
  /* drmaa_set_attribute(job, DRMAA_OUTPUT_PATH, output_path, NULL, 0); */

  /* join output/error file */
  drmaa_set_attribute(job, DRMAA_JOIN_FILES, "n", NULL, 0);

  /* transfer files */
  /* drmaa_set_attribute(job, DRMAA_TRANSFER_FILES, "ieo", NULL, 0); */
  
  /* some native SGE commands necessary */
  sprintf(attr, "-cwd -b n -q %s.q", expname);
  drmaa_set_attribute(job, DRMAA_NATIVE_SPECIFICATION, attr, NULL, 0);  

  /* print out job attributes */
  drmaa_get_attribute_names (&job_attributes, error, DRMAA_ERROR_STRING_BUFFER);

  if ( cdoVerbose )
    while ((drmaa_errno = drmaa_get_next_attr_name(job_attributes, name, DRMAA_ATTR_BUFFER)) == DRMAA_ERRNO_SUCCESS) {
      drmaa_get_attribute (job, name, value, DRMAA_ATTR_BUFFER, error,  DRMAA_ERROR_STRING_BUFFER);

      fprintf (stderr, "name: %-25s \t %s\n", name, value);
    }

  free(dir);

  return job;
}
#endif


#if defined(HAVE_LIBDRMAA)
static int drmaa_submit(const char *expname, const char *jobfilename, const char *jobname, const char *tmppath, const char *ftppath)
{
  char status[DRMAA_ERROR_STRING_BUFFER];
  char jobid[DRMAA_JOBNAME_BUFFER], jobout[DRMAA_JOBNAME_BUFFER];
  int drmaa_errno, stat;
  drmaa_job_template_t *job;
  int aborted, exited, signaled, exit_status;
  drmaa_attr_values_t *rusage = NULL;
  char usage[DRMAA_ERROR_STRING_BUFFER];
  int stdout_is_tty = 0;
  int errnum;

  { /* check character device on stdout */
    struct stat statbuf;
    fstat(1, &statbuf);
    if ( S_ISCHR(statbuf.st_mode) ) stdout_is_tty = 1;  
  }

  if ( drmaa_init(NULL, status, sizeof(status)-1) != DRMAA_ERRNO_SUCCESS )
    {
      fprintf(stderr, "drmaa_init() failed: %s\n", status);
      return 1;
    }

  /* submit some sequential jobs */

  if ( !(job = create_job_template(expname, jobfilename, jobname, tmppath)) )
    {
      fprintf(stderr, "create_job_template() failed\n");
      return 1;
    }

  while ( (drmaa_errno = drmaa_run_job(jobid, sizeof(jobid)-1, job, status, sizeof(status)-1))
	  == DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE )
    {
      fprintf(stderr, "drmaa_run_job() failed - retry: %s\n", status);
      sleep(1);
    }

  if ( drmaa_errno != DRMAA_ERRNO_SUCCESS )
    {
      fprintf(stderr, "drmaa_run_job() failed: %s\n", status);
      return 1;
    }

  if ( stdout_is_tty )
    {
      fprintf(stdout, "%s job %s ", expname, jobid);
      fprintf(stdout, "submitted ");
      fflush(stdout);
    }

  if ( stdout_is_tty )
    {
      int iwait;
      const char waitc[] = "|/-\\";

      iwait = 0;
      while ( 1 )
	{
	  sleep (1);
	
	  errnum = drmaa_job_ps(jobid, &stat, status, DRMAA_ERROR_STRING_BUFFER);
         
	  if ( errnum != DRMAA_ERRNO_SUCCESS ) break;

	  if ( stat == DRMAA_PS_QUEUED_ACTIVE ||
	       stat == DRMAA_PS_SYSTEM_ON_HOLD ||
	       stat == DRMAA_PS_USER_ON_HOLD ||
	       stat == DRMAA_PS_USER_SYSTEM_ON_HOLD )
	    {
	      fprintf(stdout, "\b\b\b\b\b\b\b\b\b\bqueued   ");
	      fprintf(stdout, "%c", (int)waitc[iwait%4]);
	      fflush(stdout);
              iwait++;
	    }
	  else
	    break;
	}

      iwait = 0;
      while ( 1 )
	{
	  sleep (1);
	
	  errnum = drmaa_job_ps(jobid, &stat, status, DRMAA_ERROR_STRING_BUFFER);

	  if ( errnum != DRMAA_ERRNO_SUCCESS ) break;

	  if ( stat == DRMAA_PS_RUNNING )
	    {
	      fprintf(stdout, "\b\b\b\b\b\b\b\b\b\brunning  ");
	      fprintf(stdout, "%c", (int)waitc[iwait%4]);
	      fflush(stdout);
              iwait++;
	    }
	  else
	    break;
	}
    }

  drmaa_delete_job_template(job, NULL, 0);

  /* wait for job */

  drmaa_errno = drmaa_wait(jobid, jobout, sizeof(jobout)-1, 
			   &stat, DRMAA_TIMEOUT_WAIT_FOREVER, &rusage, status, sizeof(status)-1);

  if ( drmaa_errno != DRMAA_ERRNO_SUCCESS )
    {
      fprintf(stderr, "drmaa_wait(%s) failed: %s\n", jobout, status);
      return 1;
    }
  
  /*
   * report how job finished 
   */
  drmaa_wifaborted(&aborted, stat, NULL, 0);
  if ( aborted )
    {
      fprintf(stderr, "job %s never ran\n", jobid);
      return 1;
    }
  else
    {
      drmaa_wifexited(&exited, stat, NULL, 0);
      if ( exited )
	{
	  drmaa_wexitstatus(&exit_status, stat, NULL, 0);
	  if ( stdout_is_tty )
	    {
	      fprintf(stdout, "\b\b\b\b\b\b\b\b\b\bfinished  \n");
	    } 

	  if ( exit_status )
	    fprintf(stdout, "%s job %s exit status %d\n", expname, jobid, exit_status);
	}
      else 
	{
	  drmaa_wifsignaled(&signaled, stat, NULL, 0);
	  if ( signaled )
	    {
	      char termsig[DRMAA_SIGNAL_BUFFER+1];
	      drmaa_wtermsig(termsig, DRMAA_SIGNAL_BUFFER, stat, NULL, 0);
	      fprintf(stderr, "job %s finished due to signal %s\n", jobid, termsig);
	    }
	  else
	    fprintf(stderr, "job %s finished with unclear conditions\n", jobid);
	}
    }

  if ( stdout_is_tty ) fprintf(stdout, "\n");

  if ( cdoVerbose )
    {
      fprintf(stderr, "Job usage:\n");
                
      while (drmaa_get_next_attr_value (rusage, usage, DRMAA_ERROR_STRING_BUFFER) == DRMAA_ERRNO_SUCCESS) {
	fprintf(stderr, "  %s\n", usage);
      }
    }
                
  drmaa_release_attr_values (rusage);

  if ( drmaa_exit(status, sizeof(status)-1) != DRMAA_ERRNO_SUCCESS )
    {
      fprintf(stderr, "drmaa_exit() failed: %s\n", status);
      return 1;
    }

  {
    char commandline[1024];
    char ftp_url[4096];
    char outname[1024];
    char errname[1024];
    int status;

    sprintf(ftp_url, "ftp://%s.zmaw.de", expname);

    sprintf(outname, "%s.o%s", jobname, jobid);
    sprintf(errname, "%s.e%s", jobname, jobid);

    status = ftpget(0, ftp_url, ftppath, outname, outname);
    if ( status == 0 )
      {
	sprintf(commandline, "cat %s | grep -v tty  | grep -v shell | grep -v SunOS | grep -v logout\n", outname);
	status = system(commandline);
      }

    status = ftpget(0, ftp_url, ftppath, errname, errname);
    if ( status == 0 )
      {
	sprintf(commandline, "cat %s | grep -v cannot | grep -v resize | grep -v rm\n", errname);
	status = system(commandline);
      }

    sprintf(commandline, "rm -f %s %s\n", outname, errname);
    status = system(commandline);
  }

  return 0;
}
#endif


int job_submit(const char *expname, const char *jobfilename, const char *jobname, const char *tmppath, const char *ftppath)
{
  int status = 0;
#if defined(HAVE_LIBDRMAA)

  status = drmaa_submit(expname, jobfilename, jobname, tmppath, ftppath);
#else
  fprintf(stderr, "DRMAA support not compiled in!\n");
#endif
  return (status);
}



#if defined(HAVE_LIBCURL)
#include <curl/curl.h>
#endif

struct FtpFile {
  const char *filename;
  FILE *stream;
};


size_t my_fwrite(void *buffer, size_t size, size_t nmemb, void *stream)
{
  struct FtpFile *out=(struct FtpFile *)stream;
  if(out && !out->stream) {
    out->stream=fopen(out->filename, "wb");
    if(!out->stream)
      return -1;
  }
  return fwrite(buffer, size, nmemb, out->stream);
}


int my_progress_func(void *stdout_is_tty,
			double t, /* dltotal */
			double d, /* dlnow */
			double ultotal,
			double ulnow)
{
  if ( *(char*)stdout_is_tty )
    {
      fprintf(stdout, "\b\b\b\b\b%4d%%", (int) (d*100/t));
      fflush(stdout);
    }

  return 0;
}


int ftpget(int flag, const char *url, const char *path, const char *target, const char *source)
{
  int status = 0;
#if defined(HAVE_LIBCURL)
  CURL *curl;
  CURLcode res;
  struct curl_slist* commands = NULL ;
  struct FtpFile ftpfile={
    NULL, /* name to store the file as if succesful */
    NULL
  };
  char filename[8192];
  char ftpcommand[1024];
  char errorbuffer[CURL_ERROR_SIZE];
  int stdout_is_tty = 0;
  char prompt[1024];

  { /* check character device on stdout */
    struct stat statbuf;
    fstat(1, &statbuf);
    if ( S_ISCHR(statbuf.st_mode) ) stdout_is_tty = 1;  
  }

  sprintf(filename, "%s%s", path, source);

  sprintf(ftpcommand, "PWD");
  commands = curl_slist_append(commands, ftpcommand) ;

  sprintf(ftpcommand, "DELE %s\n", filename);
  commands = curl_slist_append(commands, ftpcommand) ;

  sprintf(filename, "%s%s%s", url, path, source);

  if ( flag )
    {
      /*
      sprintf(prompt, "Download %-40s ", filename);
      */
      sprintf(prompt, "Download %-30s ", source);
      fprintf(stdout, "%s     ", prompt);
    }

  ftpfile.filename = target;

  curl_global_init(CURL_GLOBAL_DEFAULT);

  curl = curl_easy_init();

  if ( curl )
    {
      curl_easy_setopt(curl, CURLOPT_NETRC, CURL_NETRC_REQUIRED);

      if ( cdoVerbose )
	curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
      else
	curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);

      if ( flag )
	{
	  curl_easy_setopt(curl, CURLOPT_NOPROGRESS, FALSE);
	  curl_easy_setopt(curl, CURLOPT_PROGRESSFUNCTION, my_progress_func);
	  curl_easy_setopt(curl, CURLOPT_PROGRESSDATA, &stdout_is_tty);
	}

      curl_easy_setopt(curl, CURLOPT_FTP_SSL, CURLFTPSSL_CONTROL); 
      curl_easy_setopt(curl, CURLOPT_FTPSSLAUTH, CURLFTPAUTH_TLS);

      curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
      curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L);

      curl_easy_setopt(curl, CURLOPT_SSLKEYPASSWD, "");

      curl_easy_setopt(curl, CURLOPT_URL, filename);

      /* define callback */

      curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, my_fwrite);

      /* set a pointer to struct to pass to the callback */

      curl_easy_setopt(curl, CURLOPT_WRITEDATA, &ftpfile);

      curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, errorbuffer);

      curl_easy_setopt(curl, CURLOPT_POSTQUOTE, commands);

      res = curl_easy_perform(curl);

      curl_slist_free_all(commands);

      curl_easy_cleanup(curl);

      if ( CURLE_OK != res )
	{
	  if ( flag )
	    {
	      /* fprintf(stderr, "curl told us %d\n", res); */
	      fprintf(stderr, "%s\n", errorbuffer);
	    }
	  status = -2;
	}
    }
  else
    {
      status = -1;
    }

  if ( ftpfile.stream ) fclose(ftpfile.stream);

  curl_global_cleanup();

  if ( flag ) fprintf(stdout, "\n");
#else
  fprintf(stderr, "CURL support not compiled in!\n");
#endif

  return (status);
}


int ftprmd(const char *url, const char *path)
{
  int status = 0;
#if defined(HAVE_LIBCURL)
  CURL *curl;
  CURLcode res;
  struct curl_slist* commands = NULL ;
  char filename[8192];
  char ftpcommand[1024];
  char errorbuffer[CURL_ERROR_SIZE];

  sprintf(ftpcommand, "RMD %s\n", path);
  commands = curl_slist_append(commands, ftpcommand) ;

  /* sprintf(filename, "%s%s/tmp", url, path); */
  sprintf(filename, "%s/tmp", url); /* dummy parameter, not used! */
  /* sprintf(filename, "%s", url);*/

  curl_global_init(CURL_GLOBAL_DEFAULT);

  curl = curl_easy_init();

  if ( curl )
    {
      curl_easy_setopt(curl, CURLOPT_NETRC, CURL_NETRC_REQUIRED);

      if ( cdoVerbose )
	curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);
      else
	curl_easy_setopt(curl, CURLOPT_VERBOSE, 0L);

      curl_easy_setopt(curl, CURLOPT_FTP_SSL, CURLFTPSSL_CONTROL); 
      curl_easy_setopt(curl, CURLOPT_FTPSSLAUTH, CURLFTPAUTH_TLS);

      curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L);
      curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L);

      curl_easy_setopt(curl, CURLOPT_SSLKEYPASSWD, "");

      curl_easy_setopt(curl, CURLOPT_URL, filename);

      curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, errorbuffer);

      curl_easy_setopt(curl, CURLOPT_PREQUOTE, commands);

      res = curl_easy_perform(curl);

      curl_slist_free_all(commands);

      curl_easy_cleanup(curl);

      if ( CURLE_OK != res )
	{
	  status = -2;
	}
    }
  else
    {
      status = -1;
    }

  curl_global_cleanup();

#else
  fprintf(stderr, "CURL support not compiled in!\n");
#endif

  return (status);
}


#define  DEFAULT_CDO_REMOTE_PATH  "/client/bin/cdo"
const char *cdojobfiles  = "ftp_files";

void exp_run(int argc, char *argv[], char *cdoExpName)
{
  char commandline[65536];
  int i;
  int status;
  char jobname[1024];
  char jobfilename[1024];
  char ftp_url[4096];
  char ftpfile[1024];
  char ftppath[4096];
  char tmppath[4096];
  char ftppath0[4096];
  char tmppath0[4096];
  char tmpdir[1024];
  FILE *jobfilep, *ftpfilep;
  char *envstr;
  size_t len;
  char host[1024];

#if defined(HAVE_GETHOSTNAME)
  gethostname(host, sizeof(host));
#else
  fprintf(stderr, "Function gethostname not available!\n");
  exit(EXIT_FAILURE);
#endif

  sprintf(tmpdir, "cdo_%s_%d", host, (int) getpid());
  /*
  printf("tmpdir: >%s<\n", tmpdir);
  */

  envstr = getenv("CDO_REMOTE_PATH");
  if ( envstr )
    {
      if ( cdoVerbose )
	fprintf(stderr, "CDO_REMOTE_PATH        = %s\n", envstr);

      strcpy(commandline, envstr);
    }
  else
    {
      strcpy(commandline, DEFAULT_CDO_REMOTE_PATH);
    }

  envstr = getenv("CDO_REMOTE_TMP");
  if ( envstr )
    {
      if ( cdoVerbose )
	fprintf(stderr, "CDO_REMOTE_TMP         = %s\n", envstr);

      strcpy(tmppath, envstr);
    }
  else
    {
      sprintf(tmppath, "/%s/tmp", cdoExpName);
    }

  envstr = getenv("CDO_REMOTE_FTP");
  if ( envstr )
    {
      if ( cdoVerbose )
	fprintf(stderr, "CDO_REMOTE_FTP         = %s\n", envstr);

      strcpy(ftppath, envstr);
    }
  else
    {
      strcpy(ftppath, tmppath);
    }

  len = strlen("scratch");
  if ( strlen(tmppath) > len+1 )
    if ( memcmp(tmppath+1, "scratch", len) == 0 )
      {
	strcpy(ftppath, tmppath+len+1);
      }

  len = strlen(cdoExpName);
  if ( strlen(tmppath) > len+1 )
    if ( memcmp(tmppath+1, cdoExpName, len) == 0 )
      {
	strcpy(ftppath, tmppath+len+1);
      }

  strcat(tmppath, "/");
  strcat(ftppath, "/");

  strcpy(tmppath0, tmppath);
  strcpy(ftppath0, ftppath);

  strcat(tmppath, tmpdir);
  strcat(ftppath, tmpdir);

  strcat(tmppath, "/");
  strcat(ftppath, "/");

  if ( cdoVerbose )
    {
      fprintf(stdout, "tmppath: >%s<\n", tmppath);
      fprintf(stdout, "ftppath: >%s<\n", ftppath);
    }

  for ( i = 1; i < argc; i++ )
    {
      strcat(commandline, " ");
      strcat(commandline, argv[i]);
    }
  /*
  printf("command: >%s<\n", commandline);
  */
  sprintf(jobfilename, "%s/cdojob_%s_%d.sh", GRID_TMPDIR, host, (int) getpid());

  jobfilep = fopen(jobfilename, "w");

  if ( jobfilep == NULL )
    {
      fprintf(stderr, "Open failed on %s\n", jobfilename);
      perror(jobfilename);
      exit(EXIT_FAILURE);
    }

  fprintf(jobfilep, "#!/bin/csh\n"); /* not used !!! */
  fprintf(jobfilep, "#uname -s\n");
  fprintf(jobfilep, "#pwd\n");
  fprintf(jobfilep, "#env\n");
  fprintf(jobfilep, "#echo\n");
  fprintf(jobfilep, "#echo $SHELL\n");
  fprintf(jobfilep, "#ls -l %s\n", tmppath);
  fprintf(jobfilep, "mkdir %s\n", tmppath);
  fprintf(jobfilep, "cd %s\n", tmppath);
  fprintf(jobfilep, "#echo $LD_LIBRARY_PATH\n");
  fprintf(jobfilep, "#setenv LD_LIBRARY_PATH /opt/gridware/sge/lib/${SGE_ARCH}:$LD_LIBRARY_PATH\n");
  fprintf(jobfilep, "%s\n", commandline);
  
  fclose(jobfilep);
      
  /*sprintf(jobname, "cdo_%s", cdoExpName); */
  sprintf(jobname, "cdo_%s_%s", host, cdoExpName);

  if ( cdoVerbose)
    {
      sprintf(commandline, "cat %s\n", jobfilename);
      status = system(commandline);
    }

  status = job_submit(cdoExpName, jobfilename, jobname, tmppath0, ftppath0);
  if ( status != 0 )
    {
      fprintf(stderr, "Abort: %s job failed!\n", cdoExpName);
      exit(EXIT_FAILURE);
    }

  sprintf(commandline, "rm -f %s\n", jobfilename);
  status = system(commandline);

  sprintf(commandline, "rm -f %s\n", cdojobfiles);
  status = system(commandline);

  sprintf(ftp_url, "ftp://%s.zmaw.de", cdoExpName);

  status = ftpget(0, ftp_url, ftppath, cdojobfiles, cdojobfiles);

  if ( status == 0 )
    {
      ftpfilep = fopen(cdojobfiles, "r");
      if ( ftpfilep )
	{
	  while ( fscanf(ftpfilep, "%s\n", ftpfile) == 1 )
	    {
	      ftpget(1, ftp_url, ftppath, ftpfile, ftpfile);
	    }

	  fclose(ftpfilep);
	}
    }

  sprintf(commandline, "rm -f %s\n", cdojobfiles);
  status = system(commandline);

  status = ftprmd(ftp_url, ftppath);
}
