! This file was automatically generated, don't edit!
!
! Fortran interface for CDI library version 1.7.0
!
! Author:
! -------
! Uwe Schulzweida, MPI-MET, Hamburg,   June 2015
!

      INTEGER    CDI_MAX_NAME          
      PARAMETER (CDI_MAX_NAME           = 256)
      INTEGER    CDI_UNDEFID           
      PARAMETER (CDI_UNDEFID            = -1)
      INTEGER    CDI_GLOBAL            
      PARAMETER (CDI_GLOBAL             = -1)
!
!  Byte order
!
      INTEGER    CDI_BIGENDIAN         
      PARAMETER (CDI_BIGENDIAN          =  0)
      INTEGER    CDI_LITTLEENDIAN      
      PARAMETER (CDI_LITTLEENDIAN       =  1)
      INTEGER    CDI_REAL              
      PARAMETER (CDI_REAL               =  1)
      INTEGER    CDI_COMP              
      PARAMETER (CDI_COMP               =  2)
      INTEGER    CDI_BOTH              
      PARAMETER (CDI_BOTH               =  3)
!
!  Error identifier
!
      INTEGER    CDI_NOERR             
      PARAMETER (CDI_NOERR              =  0)
      INTEGER    CDI_EEOF              
      PARAMETER (CDI_EEOF               = -1)
      INTEGER    CDI_ESYSTEM           
      PARAMETER (CDI_ESYSTEM            = -10)
      INTEGER    CDI_EINVAL            
      PARAMETER (CDI_EINVAL             = -20)
      INTEGER    CDI_EUFTYPE           
      PARAMETER (CDI_EUFTYPE            = -21)
      INTEGER    CDI_ELIBNAVAIL        
      PARAMETER (CDI_ELIBNAVAIL         = -22)
      INTEGER    CDI_EUFSTRUCT         
      PARAMETER (CDI_EUFSTRUCT          = -23)
      INTEGER    CDI_EUNC4             
      PARAMETER (CDI_EUNC4              = -24)
      INTEGER    CDI_ELIMIT            
      PARAMETER (CDI_ELIMIT             = -99)
!
!  File types
!
      INTEGER    FILETYPE_UNDEF        
      PARAMETER (FILETYPE_UNDEF         = -1)
      INTEGER    FILETYPE_GRB          
      PARAMETER (FILETYPE_GRB           =  1)
      INTEGER    FILETYPE_GRB2         
      PARAMETER (FILETYPE_GRB2          =  2)
      INTEGER    FILETYPE_NC           
      PARAMETER (FILETYPE_NC            =  3)
      INTEGER    FILETYPE_NC2          
      PARAMETER (FILETYPE_NC2           =  4)
      INTEGER    FILETYPE_NC4          
      PARAMETER (FILETYPE_NC4           =  5)
      INTEGER    FILETYPE_NC4C         
      PARAMETER (FILETYPE_NC4C          =  6)
      INTEGER    FILETYPE_SRV          
      PARAMETER (FILETYPE_SRV           =  7)
      INTEGER    FILETYPE_EXT          
      PARAMETER (FILETYPE_EXT           =  8)
      INTEGER    FILETYPE_IEG          
      PARAMETER (FILETYPE_IEG           =  9)
!
!  Compress types
!
      INTEGER    COMPRESS_NONE         
      PARAMETER (COMPRESS_NONE          =  0)
      INTEGER    COMPRESS_SZIP         
      PARAMETER (COMPRESS_SZIP          =  1)
      INTEGER    COMPRESS_GZIP         
      PARAMETER (COMPRESS_GZIP          =  2)
      INTEGER    COMPRESS_BZIP2        
      PARAMETER (COMPRESS_BZIP2         =  3)
      INTEGER    COMPRESS_ZIP          
      PARAMETER (COMPRESS_ZIP           =  4)
      INTEGER    COMPRESS_JPEG         
      PARAMETER (COMPRESS_JPEG          =  5)
!
!  external data types
!
      INTEGER    DATATYPE_PACK         
      PARAMETER (DATATYPE_PACK          =  0)
      INTEGER    DATATYPE_PACK1        
      PARAMETER (DATATYPE_PACK1         =  1)
      INTEGER    DATATYPE_PACK2        
      PARAMETER (DATATYPE_PACK2         =  2)
      INTEGER    DATATYPE_PACK3        
      PARAMETER (DATATYPE_PACK3         =  3)
      INTEGER    DATATYPE_PACK4        
      PARAMETER (DATATYPE_PACK4         =  4)
      INTEGER    DATATYPE_PACK5        
      PARAMETER (DATATYPE_PACK5         =  5)
      INTEGER    DATATYPE_PACK6        
      PARAMETER (DATATYPE_PACK6         =  6)
      INTEGER    DATATYPE_PACK7        
      PARAMETER (DATATYPE_PACK7         =  7)
      INTEGER    DATATYPE_PACK8        
      PARAMETER (DATATYPE_PACK8         =  8)
      INTEGER    DATATYPE_PACK9        
      PARAMETER (DATATYPE_PACK9         =  9)
      INTEGER    DATATYPE_PACK10       
      PARAMETER (DATATYPE_PACK10        = 10)
      INTEGER    DATATYPE_PACK11       
      PARAMETER (DATATYPE_PACK11        = 11)
      INTEGER    DATATYPE_PACK12       
      PARAMETER (DATATYPE_PACK12        = 12)
      INTEGER    DATATYPE_PACK13       
      PARAMETER (DATATYPE_PACK13        = 13)
      INTEGER    DATATYPE_PACK14       
      PARAMETER (DATATYPE_PACK14        = 14)
      INTEGER    DATATYPE_PACK15       
      PARAMETER (DATATYPE_PACK15        = 15)
      INTEGER    DATATYPE_PACK16       
      PARAMETER (DATATYPE_PACK16        = 16)
      INTEGER    DATATYPE_PACK17       
      PARAMETER (DATATYPE_PACK17        = 17)
      INTEGER    DATATYPE_PACK18       
      PARAMETER (DATATYPE_PACK18        = 18)
      INTEGER    DATATYPE_PACK19       
      PARAMETER (DATATYPE_PACK19        = 19)
      INTEGER    DATATYPE_PACK20       
      PARAMETER (DATATYPE_PACK20        = 20)
      INTEGER    DATATYPE_PACK21       
      PARAMETER (DATATYPE_PACK21        = 21)
      INTEGER    DATATYPE_PACK22       
      PARAMETER (DATATYPE_PACK22        = 22)
      INTEGER    DATATYPE_PACK23       
      PARAMETER (DATATYPE_PACK23        = 23)
      INTEGER    DATATYPE_PACK24       
      PARAMETER (DATATYPE_PACK24        = 24)
      INTEGER    DATATYPE_PACK25       
      PARAMETER (DATATYPE_PACK25        = 25)
      INTEGER    DATATYPE_PACK26       
      PARAMETER (DATATYPE_PACK26        = 26)
      INTEGER    DATATYPE_PACK27       
      PARAMETER (DATATYPE_PACK27        = 27)
      INTEGER    DATATYPE_PACK28       
      PARAMETER (DATATYPE_PACK28        = 28)
      INTEGER    DATATYPE_PACK29       
      PARAMETER (DATATYPE_PACK29        = 29)
      INTEGER    DATATYPE_PACK30       
      PARAMETER (DATATYPE_PACK30        = 30)
      INTEGER    DATATYPE_PACK31       
      PARAMETER (DATATYPE_PACK31        = 31)
      INTEGER    DATATYPE_PACK32       
      PARAMETER (DATATYPE_PACK32        = 32)
      INTEGER    DATATYPE_CPX32        
      PARAMETER (DATATYPE_CPX32         = 64)
      INTEGER    DATATYPE_CPX64        
      PARAMETER (DATATYPE_CPX64         = 128)
      INTEGER    DATATYPE_FLT32        
      PARAMETER (DATATYPE_FLT32         = 132)
      INTEGER    DATATYPE_FLT64        
      PARAMETER (DATATYPE_FLT64         = 164)
      INTEGER    DATATYPE_INT8         
      PARAMETER (DATATYPE_INT8          = 208)
      INTEGER    DATATYPE_INT16        
      PARAMETER (DATATYPE_INT16         = 216)
      INTEGER    DATATYPE_INT32        
      PARAMETER (DATATYPE_INT32         = 232)
      INTEGER    DATATYPE_UINT8        
      PARAMETER (DATATYPE_UINT8         = 308)
      INTEGER    DATATYPE_UINT16       
      PARAMETER (DATATYPE_UINT16        = 316)
      INTEGER    DATATYPE_UINT32       
      PARAMETER (DATATYPE_UINT32        = 332)
!
!  internal data types
!
      INTEGER    DATATYPE_INT          
      PARAMETER (DATATYPE_INT           = 251)
      INTEGER    DATATYPE_FLT          
      PARAMETER (DATATYPE_FLT           = 252)
      INTEGER    DATATYPE_TXT          
      PARAMETER (DATATYPE_TXT           = 253)
      INTEGER    DATATYPE_CPX          
      PARAMETER (DATATYPE_CPX           = 254)
      INTEGER    DATATYPE_UCHAR        
      PARAMETER (DATATYPE_UCHAR         = 255)
      INTEGER    DATATYPE_LONG         
      PARAMETER (DATATYPE_LONG          = 256)
!
!  Chunks
!
      INTEGER    CHUNK_AUTO            
      PARAMETER (CHUNK_AUTO             =  1)
      INTEGER    CHUNK_GRID            
      PARAMETER (CHUNK_GRID             =  2)
      INTEGER    CHUNK_LINES           
      PARAMETER (CHUNK_LINES            =  3)
!
!  GRID types
!
      INTEGER    GRID_GENERIC          
      PARAMETER (GRID_GENERIC           =  1)
      INTEGER    GRID_GAUSSIAN         
      PARAMETER (GRID_GAUSSIAN          =  2)
      INTEGER    GRID_GAUSSIAN_REDUCED 
      PARAMETER (GRID_GAUSSIAN_REDUCED  =  3)
      INTEGER    GRID_LONLAT           
      PARAMETER (GRID_LONLAT            =  4)
      INTEGER    GRID_SPECTRAL         
      PARAMETER (GRID_SPECTRAL          =  5)
      INTEGER    GRID_FOURIER          
      PARAMETER (GRID_FOURIER           =  6)
      INTEGER    GRID_GME              
      PARAMETER (GRID_GME               =  7)
      INTEGER    GRID_TRAJECTORY       
      PARAMETER (GRID_TRAJECTORY        =  8)
      INTEGER    GRID_UNSTRUCTURED     
      PARAMETER (GRID_UNSTRUCTURED      =  9)
      INTEGER    GRID_CURVILINEAR      
      PARAMETER (GRID_CURVILINEAR       = 10)
      INTEGER    GRID_LCC              
      PARAMETER (GRID_LCC               = 11)
      INTEGER    GRID_LCC2             
      PARAMETER (GRID_LCC2              = 12)
      INTEGER    GRID_LAEA             
      PARAMETER (GRID_LAEA              = 13)
      INTEGER    GRID_SINUSOIDAL       
      PARAMETER (GRID_SINUSOIDAL        = 14)
      INTEGER    GRID_PROJECTION       
      PARAMETER (GRID_PROJECTION        = 15)
!
!  ZAXIS types
!
      INTEGER    ZAXIS_SURFACE         
      PARAMETER (ZAXIS_SURFACE          =  0)
      INTEGER    ZAXIS_GENERIC         
      PARAMETER (ZAXIS_GENERIC          =  1)
      INTEGER    ZAXIS_HYBRID          
      PARAMETER (ZAXIS_HYBRID           =  2)
      INTEGER    ZAXIS_HYBRID_HALF     
      PARAMETER (ZAXIS_HYBRID_HALF      =  3)
      INTEGER    ZAXIS_PRESSURE        
      PARAMETER (ZAXIS_PRESSURE         =  4)
      INTEGER    ZAXIS_HEIGHT          
      PARAMETER (ZAXIS_HEIGHT           =  5)
      INTEGER    ZAXIS_DEPTH_BELOW_SEA 
      PARAMETER (ZAXIS_DEPTH_BELOW_SEA  =  6)
      INTEGER    ZAXIS_DEPTH_BELOW_LAND
      PARAMETER (ZAXIS_DEPTH_BELOW_LAND =  7)
      INTEGER    ZAXIS_ISENTROPIC      
      PARAMETER (ZAXIS_ISENTROPIC       =  8)
      INTEGER    ZAXIS_TRAJECTORY      
      PARAMETER (ZAXIS_TRAJECTORY       =  9)
      INTEGER    ZAXIS_ALTITUDE        
      PARAMETER (ZAXIS_ALTITUDE         = 10)
      INTEGER    ZAXIS_SIGMA           
      PARAMETER (ZAXIS_SIGMA            = 11)
      INTEGER    ZAXIS_MEANSEA         
      PARAMETER (ZAXIS_MEANSEA          = 12)
      INTEGER    ZAXIS_TOA             
      PARAMETER (ZAXIS_TOA              = 13)
      INTEGER    ZAXIS_SEA_BOTTOM      
      PARAMETER (ZAXIS_SEA_BOTTOM       = 14)
      INTEGER    ZAXIS_ATMOSPHERE      
      PARAMETER (ZAXIS_ATMOSPHERE       = 15)
      INTEGER    ZAXIS_CLOUD_BASE      
      PARAMETER (ZAXIS_CLOUD_BASE       = 16)
      INTEGER    ZAXIS_CLOUD_TOP       
      PARAMETER (ZAXIS_CLOUD_TOP        = 17)
      INTEGER    ZAXIS_ISOTHERM_ZERO   
      PARAMETER (ZAXIS_ISOTHERM_ZERO    = 18)
      INTEGER    ZAXIS_SNOW            
      PARAMETER (ZAXIS_SNOW             = 19)
      INTEGER    ZAXIS_LAKE_BOTTOM     
      PARAMETER (ZAXIS_LAKE_BOTTOM      = 20)
      INTEGER    ZAXIS_SEDIMENT_BOTTOM 
      PARAMETER (ZAXIS_SEDIMENT_BOTTOM  = 21)
      INTEGER    ZAXIS_SEDIMENT_BOTTOM_TA
      PARAMETER (ZAXIS_SEDIMENT_BOTTOM_TA = 22)
      INTEGER    ZAXIS_SEDIMENT_BOTTOM_TW
      PARAMETER (ZAXIS_SEDIMENT_BOTTOM_TW = 23)
      INTEGER    ZAXIS_MIX_LAYER       
      PARAMETER (ZAXIS_MIX_LAYER        = 24)
      INTEGER    ZAXIS_REFERENCE       
      PARAMETER (ZAXIS_REFERENCE        = 25)
!
!  SUBTYPE types
!
      INTEGER    MAX_KV_PAIRS_MATCH    
      PARAMETER (MAX_KV_PAIRS_MATCH     = 10)
!
!  Data structure defining a key-value search, possibly with multiple
!  key-value pairs in combination.
!  
!  Currently, only multiple pairs combined by AND are supported.
!
!
!  TIME types
!
      INTEGER    TIME_CONSTANT         
      PARAMETER (TIME_CONSTANT          =  0)
      INTEGER    TIME_VARIABLE         
      PARAMETER (TIME_VARIABLE          =  1)
!
!  TSTEP types
!
      INTEGER    TSTEP_CONSTANT        
      PARAMETER (TSTEP_CONSTANT         =  0)
      INTEGER    TSTEP_INSTANT         
      PARAMETER (TSTEP_INSTANT          =  1)
      INTEGER    TSTEP_AVG             
      PARAMETER (TSTEP_AVG              =  2)
      INTEGER    TSTEP_ACCUM           
      PARAMETER (TSTEP_ACCUM            =  3)
      INTEGER    TSTEP_MAX             
      PARAMETER (TSTEP_MAX              =  4)
      INTEGER    TSTEP_MIN             
      PARAMETER (TSTEP_MIN              =  5)
      INTEGER    TSTEP_DIFF            
      PARAMETER (TSTEP_DIFF             =  6)
      INTEGER    TSTEP_RMS             
      PARAMETER (TSTEP_RMS              =  7)
      INTEGER    TSTEP_SD              
      PARAMETER (TSTEP_SD               =  8)
      INTEGER    TSTEP_COV             
      PARAMETER (TSTEP_COV              =  9)
      INTEGER    TSTEP_RATIO           
      PARAMETER (TSTEP_RATIO            = 10)
      INTEGER    TSTEP_RANGE           
      PARAMETER (TSTEP_RANGE            = 11)
      INTEGER    TSTEP_INSTANT2        
      PARAMETER (TSTEP_INSTANT2         = 12)
      INTEGER    TSTEP_INSTANT3        
      PARAMETER (TSTEP_INSTANT3         = 13)
!
!  TAXIS types
!
      INTEGER    TAXIS_ABSOLUTE        
      PARAMETER (TAXIS_ABSOLUTE         =  1)
      INTEGER    TAXIS_RELATIVE        
      PARAMETER (TAXIS_RELATIVE         =  2)
      INTEGER    TAXIS_FORECAST        
      PARAMETER (TAXIS_FORECAST         =  3)
!
!  TUNIT types
!
      INTEGER    TUNIT_SECOND          
      PARAMETER (TUNIT_SECOND           =  1)
      INTEGER    TUNIT_MINUTE          
      PARAMETER (TUNIT_MINUTE           =  2)
      INTEGER    TUNIT_QUARTER         
      PARAMETER (TUNIT_QUARTER          =  3)
      INTEGER    TUNIT_30MINUTES       
      PARAMETER (TUNIT_30MINUTES        =  4)
      INTEGER    TUNIT_HOUR            
      PARAMETER (TUNIT_HOUR             =  5)
      INTEGER    TUNIT_3HOURS          
      PARAMETER (TUNIT_3HOURS           =  6)
      INTEGER    TUNIT_6HOURS          
      PARAMETER (TUNIT_6HOURS           =  7)
      INTEGER    TUNIT_12HOURS         
      PARAMETER (TUNIT_12HOURS          =  8)
      INTEGER    TUNIT_DAY             
      PARAMETER (TUNIT_DAY              =  9)
      INTEGER    TUNIT_MONTH           
      PARAMETER (TUNIT_MONTH            = 10)
      INTEGER    TUNIT_YEAR            
      PARAMETER (TUNIT_YEAR             = 11)
!
!  CALENDAR types
!
      INTEGER    CALENDAR_STANDARD     
      PARAMETER (CALENDAR_STANDARD      =  0)
      INTEGER    CALENDAR_PROLEPTIC    
      PARAMETER (CALENDAR_PROLEPTIC     =  1)
      INTEGER    CALENDAR_360DAYS      
      PARAMETER (CALENDAR_360DAYS       =  2)
      INTEGER    CALENDAR_365DAYS      
      PARAMETER (CALENDAR_365DAYS       =  3)
      INTEGER    CALENDAR_366DAYS      
      PARAMETER (CALENDAR_366DAYS       =  4)
      INTEGER    CALENDAR_NONE         
      PARAMETER (CALENDAR_NONE          =  5)
!
!  number of unsigned char needed to store UUID
!
      INTEGER    CDI_UUID_SIZE         
      PARAMETER (CDI_UUID_SIZE          = 16)
!
!  Structs that are used to return data to the user
!
!
!  Opaque types
!
!
!  CDI control routines
!
!                     cdiReset
      EXTERNAL        cdiReset

      CHARACTER(80)   cdiStringError
!                                    (INTEGER         cdiErrno)
      EXTERNAL        cdiStringError

!                     cdiDebug
!                                    (INTEGER         debug)
      EXTERNAL        cdiDebug

      CHARACTER(80)   cdiLibraryVersion
      EXTERNAL        cdiLibraryVersion

!                     cdiPrintVersion
      EXTERNAL        cdiPrintVersion

      INTEGER         cdiHaveFiletype
!                                    (INTEGER         filetype)
      EXTERNAL        cdiHaveFiletype

!                     cdiDefMissval
!                                    (DOUBLEPRECISION missval)
      EXTERNAL        cdiDefMissval

      DOUBLEPRECISION cdiInqMissval
      EXTERNAL        cdiInqMissval

!                     cdiDefGlobal
!                                    (CHARACTER*(*)   string,
!                                     INTEGER         val)
      EXTERNAL        cdiDefGlobal

      INTEGER         namespaceNew
      EXTERNAL        namespaceNew

!                     namespaceSetActive
!                                    (INTEGER         namespaceID)
      EXTERNAL        namespaceSetActive

!                     namespaceDelete
!                                    (INTEGER         namespaceID)
      EXTERNAL        namespaceDelete

!
!  CDI converter routines
!
!
!  parameter
!
!                     cdiParamToString
!                                    (INTEGER         param,
!                                     CHARACTER*(*)   paramstr,
!                                     INTEGER         maxlen)
      EXTERNAL        cdiParamToString

!                     cdiDecodeParam
!                                    (INTEGER         param,
!                                     INTEGER         pnum,
!                                     INTEGER         pcat,
!                                     INTEGER         pdis)
      EXTERNAL        cdiDecodeParam

      INTEGER         cdiEncodeParam
!                                    (INTEGER         pnum,
!                                     INTEGER         pcat,
!                                     INTEGER         pdis)
      EXTERNAL        cdiEncodeParam

!
!  date format:  YYYYMMDD
!
!
!  time format:    hhmmss
!
!                     cdiDecodeDate
!                                    (INTEGER         date,
!                                     INTEGER         year,
!                                     INTEGER         month,
!                                     INTEGER         day)
      EXTERNAL        cdiDecodeDate

      INTEGER         cdiEncodeDate
!                                    (INTEGER         year,
!                                     INTEGER         month,
!                                     INTEGER         day)
      EXTERNAL        cdiEncodeDate

!                     cdiDecodeTime
!                                    (INTEGER         time,
!                                     INTEGER         hour,
!                                     INTEGER         minute,
!                                     INTEGER         second)
      EXTERNAL        cdiDecodeTime

      INTEGER         cdiEncodeTime
!                                    (INTEGER         hour,
!                                     INTEGER         minute,
!                                     INTEGER         second)
      EXTERNAL        cdiEncodeTime

!
!  STREAM control routines
!
      INTEGER         cdiGetFiletype
!                                    (CHARACTER*(*)   path,
!                                     INTEGER         byteorder)
      EXTERNAL        cdiGetFiletype

      INTEGER         streamOpenRead
!                                    (CHARACTER*(*)   path)
      EXTERNAL        streamOpenRead

      INTEGER         streamOpenWrite
!                                    (CHARACTER*(*)   path,
!                                     INTEGER         filetype)
      EXTERNAL        streamOpenWrite

      INTEGER         streamOpenAppend
!                                    (CHARACTER*(*)   path)
      EXTERNAL        streamOpenAppend

!                     streamClose
!                                    (INTEGER         streamID)
      EXTERNAL        streamClose

!                     streamSync
!                                    (INTEGER         streamID)
      EXTERNAL        streamSync

!                     streamDefVlist
!                                    (INTEGER         streamID,
!                                     INTEGER         vlistID)
      EXTERNAL        streamDefVlist

      INTEGER         streamInqVlist
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqVlist

!
!  PIO:
!
      INTEGER         streamInqVlistIDorig
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqVlistIDorig

      INTEGER         streamInqFiletype
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqFiletype

!                     streamDefByteorder
!                                    (INTEGER         streamID,
!                                     INTEGER         byteorder)
      EXTERNAL        streamDefByteorder

      INTEGER         streamInqByteorder
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqByteorder

!                     streamDefCompType
!                                    (INTEGER         streamID,
!                                     INTEGER         comptype)
      EXTERNAL        streamDefCompType

      INTEGER         streamInqCompType
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqCompType

!                     streamDefCompLevel
!                                    (INTEGER         streamID,
!                                     INTEGER         complevel)
      EXTERNAL        streamDefCompLevel

      INTEGER         streamInqCompLevel
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqCompLevel

      INTEGER         streamDefTimestep
!                                    (INTEGER         streamID,
!                                     INTEGER         tsID)
      EXTERNAL        streamDefTimestep

      INTEGER         streamInqTimestep
!                                    (INTEGER         streamID,
!                                     INTEGER         tsID)
      EXTERNAL        streamInqTimestep

      INTEGER         streamInqCurTimestepID
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqCurTimestepID

      CHARACTER(80)   streamFilename
!                                    (INTEGER         streamID)
      EXTERNAL        streamFilename

      CHARACTER(80)   streamFilesuffix
!                                    (INTEGER         filetype)
      EXTERNAL        streamFilesuffix

      INTEGER         streamInqNvars
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqNvars

!
!  STREAM var I/O routines
!
!                     streamWriteVar
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     DOUBLEPRECISION data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamWriteVar

!                     streamWriteVarF
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     REAL            data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamWriteVarF

!                     streamReadVar
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     DOUBLEPRECISION data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamReadVar

!                     streamReadVarF
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     REAL            data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamReadVarF

!                     streamWriteVarSlice
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     INTEGER         levelID,
!                                     DOUBLEPRECISION data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamWriteVarSlice

!                     streamWriteVarSliceF
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     INTEGER         levelID,
!                                     REAL            data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamWriteVarSliceF

!                     streamReadVarSlice
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     INTEGER         levelID,
!                                     DOUBLEPRECISION data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamReadVarSlice

!                     streamReadVarSliceF
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     INTEGER         levelID,
!                                     REAL            data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamReadVarSliceF

!                     streamWriteVarChunk
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     INTEGER         rect,
!                                     DOUBLEPRECISION data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamWriteVarChunk

!
!  STREAM record I/O routines
!
!                     streamDefRecord
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     INTEGER         levelID)
      EXTERNAL        streamDefRecord

!                     streamInqRecord
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     INTEGER         levelID)
      EXTERNAL        streamInqRecord

!                     streamWriteRecord
!                                    (INTEGER         streamID,
!                                     DOUBLEPRECISION data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamWriteRecord

!                     streamWriteRecordF
!                                    (INTEGER         streamID,
!                                     REAL            data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamWriteRecordF

!                     streamReadRecord
!                                    (INTEGER         streamID,
!                                     DOUBLEPRECISION data_vec,
!                                     INTEGER         nmiss)
      EXTERNAL        streamReadRecord

!                     streamCopyRecord
!                                    (INTEGER         streamIDdest,
!                                     INTEGER         streamIDsrc)
      EXTERNAL        streamCopyRecord

!
!  File driven I/O (may yield better performance than using the streamXXX functions)
!
!
!  VLIST routines
!
      INTEGER         vlistCreate
      EXTERNAL        vlistCreate

!                     vlistDestroy
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistDestroy

      INTEGER         vlistDuplicate
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistDuplicate

!                     vlistCopy
!                                    (INTEGER         vlistID2,
!                                     INTEGER         vlistID1)
      EXTERNAL        vlistCopy

!                     vlistCopyFlag
!                                    (INTEGER         vlistID2,
!                                     INTEGER         vlistID1)
      EXTERNAL        vlistCopyFlag

!                     vlistClearFlag
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistClearFlag

!                     vlistCat
!                                    (INTEGER         vlistID2,
!                                     INTEGER         vlistID1)
      EXTERNAL        vlistCat

!                     vlistMerge
!                                    (INTEGER         vlistID2,
!                                     INTEGER         vlistID1)
      EXTERNAL        vlistMerge

!                     vlistPrint
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistPrint

      INTEGER         vlistNumber
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistNumber

      INTEGER         vlistNvars
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistNvars

      INTEGER         vlistNgrids
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistNgrids

      INTEGER         vlistNzaxis
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistNzaxis

      INTEGER         vlistNsubtypes
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistNsubtypes

!                     vlistDefNtsteps
!                                    (INTEGER         vlistID,
!                                     INTEGER         nts)
      EXTERNAL        vlistDefNtsteps

      INTEGER         vlistNtsteps
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistNtsteps

      INTEGER         vlistGridsizeMax
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistGridsizeMax

      INTEGER         vlistGrid
!                                    (INTEGER         vlistID,
!                                     INTEGER         index)
      EXTERNAL        vlistGrid

      INTEGER         vlistGridIndex
!                                    (INTEGER         vlistID,
!                                     INTEGER         gridID)
      EXTERNAL        vlistGridIndex

!                     vlistChangeGridIndex
!                                    (INTEGER         vlistID,
!                                     INTEGER         index,
!                                     INTEGER         gridID)
      EXTERNAL        vlistChangeGridIndex

!                     vlistChangeGrid
!                                    (INTEGER         vlistID,
!                                     INTEGER         gridID1,
!                                     INTEGER         gridID2)
      EXTERNAL        vlistChangeGrid

      INTEGER         vlistZaxis
!                                    (INTEGER         vlistID,
!                                     INTEGER         index)
      EXTERNAL        vlistZaxis

      INTEGER         vlistZaxisIndex
!                                    (INTEGER         vlistID,
!                                     INTEGER         zaxisID)
      EXTERNAL        vlistZaxisIndex

!                     vlistChangeZaxisIndex
!                                    (INTEGER         vlistID,
!                                     INTEGER         index,
!                                     INTEGER         zaxisID)
      EXTERNAL        vlistChangeZaxisIndex

!                     vlistChangeZaxis
!                                    (INTEGER         vlistID,
!                                     INTEGER         zaxisID1,
!                                     INTEGER         zaxisID2)
      EXTERNAL        vlistChangeZaxis

      INTEGER         vlistNrecs
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistNrecs

      INTEGER         vlistSubtype
!                                    (INTEGER         vlistID,
!                                     INTEGER         index)
      EXTERNAL        vlistSubtype

      INTEGER         vlistSubtypeIndex
!                                    (INTEGER         vlistID,
!                                     INTEGER         subtypeID)
      EXTERNAL        vlistSubtypeIndex

!                     vlistDefTaxis
!                                    (INTEGER         vlistID,
!                                     INTEGER         taxisID)
      EXTERNAL        vlistDefTaxis

      INTEGER         vlistInqTaxis
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistInqTaxis

!                     vlistDefTable
!                                    (INTEGER         vlistID,
!                                     INTEGER         tableID)
      EXTERNAL        vlistDefTable

      INTEGER         vlistInqTable
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistInqTable

!                     vlistDefInstitut
!                                    (INTEGER         vlistID,
!                                     INTEGER         instID)
      EXTERNAL        vlistDefInstitut

      INTEGER         vlistInqInstitut
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistInqInstitut

!                     vlistDefModel
!                                    (INTEGER         vlistID,
!                                     INTEGER         modelID)
      EXTERNAL        vlistDefModel

      INTEGER         vlistInqModel
!                                    (INTEGER         vlistID)
      EXTERNAL        vlistInqModel

!
!  VLIST VAR routines
!
      INTEGER         vlistDefVarTiles
!                                    (INTEGER         vlistID,
!                                     INTEGER         gridID,
!                                     INTEGER         zaxisID,
!                                     INTEGER         tsteptype,
!                                     INTEGER         tilesetID)
      EXTERNAL        vlistDefVarTiles

      INTEGER         vlistDefVar
!                                    (INTEGER         vlistID,
!                                     INTEGER         gridID,
!                                     INTEGER         zaxisID,
!                                     INTEGER         tsteptype)
      EXTERNAL        vlistDefVar

!                     vlistChangeVarGrid
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         gridID)
      EXTERNAL        vlistChangeVarGrid

!                     vlistChangeVarZaxis
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         zaxisID)
      EXTERNAL        vlistChangeVarZaxis

!                     vlistInqVar
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         gridID,
!                                     INTEGER         zaxisID,
!                                     INTEGER         tsteptype)
      EXTERNAL        vlistInqVar

      INTEGER         vlistInqVarGrid
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarGrid

      INTEGER         vlistInqVarZaxis
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarZaxis

!
!  used in MPIOM
!
      INTEGER         vlistInqVarID
!                                    (INTEGER         vlistID,
!                                     INTEGER         code)
      EXTERNAL        vlistInqVarID

!                     vlistDefVarTsteptype
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         tsteptype)
      EXTERNAL        vlistDefVarTsteptype

      INTEGER         vlistInqVarTsteptype
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarTsteptype

!                     vlistDefVarCompType
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         comptype)
      EXTERNAL        vlistDefVarCompType

      INTEGER         vlistInqVarCompType
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarCompType

!                     vlistDefVarCompLevel
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         complevel)
      EXTERNAL        vlistDefVarCompLevel

      INTEGER         vlistInqVarCompLevel
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarCompLevel

!                     vlistDefVarParam
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         param)
      EXTERNAL        vlistDefVarParam

      INTEGER         vlistInqVarParam
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarParam

!                     vlistDefVarCode
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         code)
      EXTERNAL        vlistDefVarCode

      INTEGER         vlistInqVarCode
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarCode

!                     vlistDefVarDatatype
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         datatype)
      EXTERNAL        vlistDefVarDatatype

      INTEGER         vlistInqVarDatatype
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarDatatype

!                     vlistDefVarChunkType
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         chunktype)
      EXTERNAL        vlistDefVarChunkType

      INTEGER         vlistInqVarChunkType
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarChunkType

!                     vlistDefVarXYZ
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         xyz)
      EXTERNAL        vlistDefVarXYZ

      INTEGER         vlistInqVarXYZ
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarXYZ

      INTEGER         vlistInqVarNumber
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarNumber

!                     vlistDefVarInstitut
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         instID)
      EXTERNAL        vlistDefVarInstitut

      INTEGER         vlistInqVarInstitut
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarInstitut

!                     vlistDefVarModel
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         modelID)
      EXTERNAL        vlistDefVarModel

      INTEGER         vlistInqVarModel
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarModel

!                     vlistDefVarTable
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         tableID)
      EXTERNAL        vlistDefVarTable

      INTEGER         vlistInqVarTable
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarTable

!                     vlistDefVarName
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        vlistDefVarName

!                     vlistInqVarName
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        vlistInqVarName

      CHARACTER(80)   vlistCopyVarName
!                                    (INTEGER         vlistId,
!                                     INTEGER         varId)
      EXTERNAL        vlistCopyVarName

!                     vlistDefVarStdname
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   stdname)
      EXTERNAL        vlistDefVarStdname

!                     vlistInqVarStdname
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   stdname)
      EXTERNAL        vlistInqVarStdname

!                     vlistDefVarLongname
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   longname)
      EXTERNAL        vlistDefVarLongname

!                     vlistInqVarLongname
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   longname)
      EXTERNAL        vlistInqVarLongname

!                     vlistDefVarUnits
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   units)
      EXTERNAL        vlistDefVarUnits

!                     vlistInqVarUnits
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   units)
      EXTERNAL        vlistInqVarUnits

!                     vlistDefVarMissval
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     DOUBLEPRECISION missval)
      EXTERNAL        vlistDefVarMissval

      DOUBLEPRECISION vlistInqVarMissval
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarMissval

!                     vlistDefVarExtra
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   extra)
      EXTERNAL        vlistDefVarExtra

!                     vlistInqVarExtra
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   extra)
      EXTERNAL        vlistInqVarExtra

!                     vlistDefVarScalefactor
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     DOUBLEPRECISION scalefactor)
      EXTERNAL        vlistDefVarScalefactor

      DOUBLEPRECISION vlistInqVarScalefactor
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarScalefactor

!                     vlistDefVarAddoffset
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     DOUBLEPRECISION addoffset)
      EXTERNAL        vlistDefVarAddoffset

      DOUBLEPRECISION vlistInqVarAddoffset
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarAddoffset

!                     vlistDefVarTimave
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         timave)
      EXTERNAL        vlistDefVarTimave

      INTEGER         vlistInqVarTimave
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarTimave

!                     vlistDefVarTimaccu
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         timaccu)
      EXTERNAL        vlistDefVarTimaccu

      INTEGER         vlistInqVarTimaccu
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarTimaccu

!                     vlistDefVarTypeOfGeneratingProcess
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         typeOfGeneratingProcess)
      EXTERNAL        vlistDefVarTypeOfGeneratingProcess

      INTEGER         vlistInqVarTypeOfGeneratingProcess
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarTypeOfGeneratingProcess

!                     vlistDefVarProductDefinitionTemplate
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         productDefinitionTemplate)
      EXTERNAL        vlistDefVarProductDefinitionTemplate

      INTEGER         vlistInqVarProductDefinitionTemplate
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarProductDefinitionTemplate

      INTEGER         vlistInqVarSize
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarSize

!                     vlistDefIndex
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         levID,
!                                     INTEGER         index)
      EXTERNAL        vlistDefIndex

      INTEGER         vlistInqIndex
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         levID)
      EXTERNAL        vlistInqIndex

!                     vlistDefFlag
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         levID,
!                                     INTEGER         flag)
      EXTERNAL        vlistDefFlag

      INTEGER         vlistInqFlag
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         levID)
      EXTERNAL        vlistInqFlag

      INTEGER         vlistFindVar
!                                    (INTEGER         vlistID,
!                                     INTEGER         fvarID)
      EXTERNAL        vlistFindVar

      INTEGER         vlistFindLevel
!                                    (INTEGER         vlistID,
!                                     INTEGER         fvarID,
!                                     INTEGER         flevelID)
      EXTERNAL        vlistFindLevel

      INTEGER         vlistMergedVar
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistMergedVar

      INTEGER         vlistMergedLevel
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         levelID)
      EXTERNAL        vlistMergedLevel

!
!  Ensemble info routines
!
!                     vlistDefVarEnsemble
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         ensID,
!                                     INTEGER         ensCount,
!                                     INTEGER         forecast_type)
      EXTERNAL        vlistDefVarEnsemble

      INTEGER         vlistInqVarEnsemble
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         ensID,
!                                     INTEGER         ensCount,
!                                     INTEGER         forecast_type)
      EXTERNAL        vlistInqVarEnsemble

!                     cdiClearAdditionalKeys
      EXTERNAL        cdiClearAdditionalKeys

!                     cdiDefAdditionalKey
!                                    (CHARACTER*(*)   string)
      EXTERNAL        cdiDefAdditionalKey

!                     vlistDefVarIntKey
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         value)
      EXTERNAL        vlistDefVarIntKey

!                     vlistDefVarDblKey
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name,
!                                     DOUBLEPRECISION value)
      EXTERNAL        vlistDefVarDblKey

      INTEGER         vlistHasVarKey
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        vlistHasVarKey

      DOUBLEPRECISION vlistInqVarDblKey
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        vlistInqVarDblKey

      INTEGER         vlistInqVarIntKey
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        vlistInqVarIntKey

!
!  VLIST attributes
!
      INTEGER         vlistInqNatts
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         nattsp)
      EXTERNAL        vlistInqNatts

      INTEGER         vlistInqAtt
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     INTEGER         attrnum,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         typep,
!                                     INTEGER         lenp)
      EXTERNAL        vlistInqAtt

      INTEGER         vlistDelAtt
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        vlistDelAtt

      INTEGER         vlistDefAttInt
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         type,
!                                     INTEGER         len,
!                                     INTEGER         ip_vec)
      EXTERNAL        vlistDefAttInt

      INTEGER         vlistDefAttFlt
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         type,
!                                     INTEGER         len,
!                                     DOUBLEPRECISION dp_vec)
      EXTERNAL        vlistDefAttFlt

      INTEGER         vlistDefAttTxt
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         len,
!                                     CHARACTER*(*)   tp_cbuf)
      EXTERNAL        vlistDefAttTxt

      INTEGER         vlistInqAttInt
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         mlen,
!                                     INTEGER         ip_vec)
      EXTERNAL        vlistInqAttInt

      INTEGER         vlistInqAttFlt
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         mlen,
!                                     DOUBLEPRECISION dp_vec)
      EXTERNAL        vlistInqAttFlt

      INTEGER         vlistInqAttTxt
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         mlen,
!                                     CHARACTER*(*)   tp_cbuf)
      EXTERNAL        vlistInqAttTxt

!
!  GRID routines
!
!                     gridName
!                                    (INTEGER         gridtype,
!                                     CHARACTER*(*)   gridname)
      EXTERNAL        gridName

      CHARACTER(80)   gridNamePtr
!                                    (INTEGER         gridtype)
      EXTERNAL        gridNamePtr

!                     gridCompress
!                                    (INTEGER         gridID)
      EXTERNAL        gridCompress

!                     gridDefMaskGME
!                                    (INTEGER         gridID,
!                                     INTEGER         mask_vec)
      EXTERNAL        gridDefMaskGME

      INTEGER         gridInqMaskGME
!                                    (INTEGER         gridID,
!                                     INTEGER         mask_vec)
      EXTERNAL        gridInqMaskGME

!                     gridDefMask
!                                    (INTEGER         gridID,
!                                     INTEGER         mask_vec)
      EXTERNAL        gridDefMask

      INTEGER         gridInqMask
!                                    (INTEGER         gridID,
!                                     INTEGER         mask_vec)
      EXTERNAL        gridInqMask

!                     gridPrint
!                                    (INTEGER         gridID,
!                                     INTEGER         index,
!                                     INTEGER         opt)
      EXTERNAL        gridPrint

      INTEGER         gridCreate
!                                    (INTEGER         gridtype,
!                                     INTEGER         size)
      EXTERNAL        gridCreate

!                     gridDestroy
!                                    (INTEGER         gridID)
      EXTERNAL        gridDestroy

      INTEGER         gridDuplicate
!                                    (INTEGER         gridID)
      EXTERNAL        gridDuplicate

      INTEGER         gridInqType
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqType

      INTEGER         gridInqSize
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqSize

!                     gridDefXsize
!                                    (INTEGER         gridID,
!                                     INTEGER         xsize)
      EXTERNAL        gridDefXsize

      INTEGER         gridInqXsize
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqXsize

!                     gridDefYsize
!                                    (INTEGER         gridID,
!                                     INTEGER         ysize)
      EXTERNAL        gridDefYsize

      INTEGER         gridInqYsize
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqYsize

!                     gridDefNP
!                                    (INTEGER         gridID,
!                                     INTEGER         np)
      EXTERNAL        gridDefNP

      INTEGER         gridInqNP
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqNP

!                     gridDefXvals
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION xvals_vec)
      EXTERNAL        gridDefXvals

      INTEGER         gridInqXvals
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION xvals_vec)
      EXTERNAL        gridInqXvals

!                     gridDefYvals
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION yvals_vec)
      EXTERNAL        gridDefYvals

      INTEGER         gridInqYvals
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION yvals_vec)
      EXTERNAL        gridInqYvals

!                     gridDefXname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   xname)
      EXTERNAL        gridDefXname

!                     gridInqXname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   xname)
      EXTERNAL        gridInqXname

!                     gridDefXlongname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   xlongname)
      EXTERNAL        gridDefXlongname

!                     gridInqXlongname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   xlongname)
      EXTERNAL        gridInqXlongname

!                     gridDefXunits
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   xunits)
      EXTERNAL        gridDefXunits

!                     gridInqXunits
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   xunits)
      EXTERNAL        gridInqXunits

!                     gridDefYname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   yname)
      EXTERNAL        gridDefYname

!                     gridInqYname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   yname)
      EXTERNAL        gridInqYname

!                     gridDefYlongname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   ylongname)
      EXTERNAL        gridDefYlongname

!                     gridInqYlongname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   ylongname)
      EXTERNAL        gridInqYlongname

!                     gridDefYunits
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   yunits)
      EXTERNAL        gridDefYunits

!                     gridInqYunits
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   yunits)
      EXTERNAL        gridInqYunits

!                     gridInqXstdname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   xstdname)
      EXTERNAL        gridInqXstdname

!                     gridInqYstdname
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   ystdname)
      EXTERNAL        gridInqYstdname

!                     gridDefPrec
!                                    (INTEGER         gridID,
!                                     INTEGER         prec)
      EXTERNAL        gridDefPrec

      INTEGER         gridInqPrec
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqPrec

      DOUBLEPRECISION gridInqXval
!                                    (INTEGER         gridID,
!                                     INTEGER         index)
      EXTERNAL        gridInqXval

      DOUBLEPRECISION gridInqYval
!                                    (INTEGER         gridID,
!                                     INTEGER         index)
      EXTERNAL        gridInqYval

      DOUBLEPRECISION gridInqXinc
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqXinc

      DOUBLEPRECISION gridInqYinc
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqYinc

      INTEGER         gridIsCircular
!                                    (INTEGER         gridID)
      EXTERNAL        gridIsCircular

      INTEGER         gridIsRotated
!                                    (INTEGER         gridID)
      EXTERNAL        gridIsRotated

!                     gridDefXpole
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION xpole)
      EXTERNAL        gridDefXpole

      DOUBLEPRECISION gridInqXpole
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqXpole

!                     gridDefYpole
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION ypole)
      EXTERNAL        gridDefYpole

      DOUBLEPRECISION gridInqYpole
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqYpole

!                     gridDefAngle
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION angle)
      EXTERNAL        gridDefAngle

      DOUBLEPRECISION gridInqAngle
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqAngle

      INTEGER         gridInqTrunc
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqTrunc

!                     gridDefTrunc
!                                    (INTEGER         gridID,
!                                     INTEGER         trunc)
      EXTERNAL        gridDefTrunc

!
!  Hexagonal GME grid
!
!                     gridDefGMEnd
!                                    (INTEGER         gridID,
!                                     INTEGER         nd)
      EXTERNAL        gridDefGMEnd

      INTEGER         gridInqGMEnd
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqGMEnd

!                     gridDefGMEni
!                                    (INTEGER         gridID,
!                                     INTEGER         ni)
      EXTERNAL        gridDefGMEni

      INTEGER         gridInqGMEni
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqGMEni

!                     gridDefGMEni2
!                                    (INTEGER         gridID,
!                                     INTEGER         ni2)
      EXTERNAL        gridDefGMEni2

      INTEGER         gridInqGMEni2
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqGMEni2

!                     gridDefGMEni3
!                                    (INTEGER         gridID,
!                                     INTEGER         ni3)
      EXTERNAL        gridDefGMEni3

      INTEGER         gridInqGMEni3
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqGMEni3

!
!  Reference of an unstructured grid
!
!                     gridDefNumber
!                                    (INTEGER         gridID,
!                                     INTEGER         number)
      EXTERNAL        gridDefNumber

      INTEGER         gridInqNumber
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqNumber

!                     gridDefPosition
!                                    (INTEGER         gridID,
!                                     INTEGER         position)
      EXTERNAL        gridDefPosition

      INTEGER         gridInqPosition
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqPosition

!                     gridDefReference
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   reference)
      EXTERNAL        gridDefReference

      INTEGER         gridInqReference
!                                    (INTEGER         gridID,
!                                     CHARACTER*(*)   reference)
      EXTERNAL        gridInqReference

!                     gridDefUUID
!                                    (INTEGER         gridID,
!                                     INTEGER*1(16)   uuid)
      EXTERNAL        gridDefUUID

!                     gridInqUUID
!                                    (INTEGER         gridID,
!                                     INTEGER*1(16)   uuid)
      EXTERNAL        gridInqUUID

!
!  Lambert Conformal Conic grid (GRIB version)
!
!                     gridDefLCC
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION originLon,
!                                     DOUBLEPRECISION originLat,
!                                     DOUBLEPRECISION lonParY,
!                                     DOUBLEPRECISION lat1,
!                                     DOUBLEPRECISION lat2,
!                                     DOUBLEPRECISION xinc,
!                                     DOUBLEPRECISION yinc,
!                                     INTEGER         projflag,
!                                     INTEGER         scanflag)
      EXTERNAL        gridDefLCC

!                     gridInqLCC
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION originLon,
!                                     DOUBLEPRECISION originLat,
!                                     DOUBLEPRECISION lonParY,
!                                     DOUBLEPRECISION lat1,
!                                     DOUBLEPRECISION lat2,
!                                     DOUBLEPRECISION xinc,
!                                     DOUBLEPRECISION yinc,
!                                     INTEGER         projflag,
!                                     INTEGER         scanflag)
      EXTERNAL        gridInqLCC

!
!  Lambert Conformal Conic 2 grid (PROJ version)
!
!                     gridDefLcc2
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION earth_radius,
!                                     DOUBLEPRECISION lon_0,
!                                     DOUBLEPRECISION lat_0,
!                                     DOUBLEPRECISION lat_1,
!                                     DOUBLEPRECISION lat_2)
      EXTERNAL        gridDefLcc2

!                     gridInqLcc2
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION earth_radius,
!                                     DOUBLEPRECISION lon_0,
!                                     DOUBLEPRECISION lat_0,
!                                     DOUBLEPRECISION lat_1,
!                                     DOUBLEPRECISION lat_2)
      EXTERNAL        gridInqLcc2

!
!  Lambert Azimuthal Equal Area grid
!
!                     gridDefLaea
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION earth_radius,
!                                     DOUBLEPRECISION lon_0,
!                                     DOUBLEPRECISION lat_0)
      EXTERNAL        gridDefLaea

!                     gridInqLaea
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION earth_radius,
!                                     DOUBLEPRECISION lon_0,
!                                     DOUBLEPRECISION lat_0)
      EXTERNAL        gridInqLaea

!                     gridDefArea
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION area_vec)
      EXTERNAL        gridDefArea

!                     gridInqArea
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION area_vec)
      EXTERNAL        gridInqArea

      INTEGER         gridHasArea
!                                    (INTEGER         gridID)
      EXTERNAL        gridHasArea

!                     gridDefNvertex
!                                    (INTEGER         gridID,
!                                     INTEGER         nvertex)
      EXTERNAL        gridDefNvertex

      INTEGER         gridInqNvertex
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqNvertex

!                     gridDefXbounds
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION xbounds_vec)
      EXTERNAL        gridDefXbounds

      INTEGER         gridInqXbounds
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION xbounds_vec)
      EXTERNAL        gridInqXbounds

!                     gridDefYbounds
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION ybounds_vec)
      EXTERNAL        gridDefYbounds

      INTEGER         gridInqYbounds
!                                    (INTEGER         gridID,
!                                     DOUBLEPRECISION ybounds_vec)
      EXTERNAL        gridInqYbounds

!                     gridDefRowlon
!                                    (INTEGER         gridID,
!                                     INTEGER         nrowlon,
!                                     INTEGER         rowlon_vec)
      EXTERNAL        gridDefRowlon

!                     gridInqRowlon
!                                    (INTEGER         gridID,
!                                     INTEGER         rowlon_vec)
      EXTERNAL        gridInqRowlon

!                     gridChangeType
!                                    (INTEGER         gridID,
!                                     INTEGER         gridtype)
      EXTERNAL        gridChangeType

!                     gridDefComplexPacking
!                                    (INTEGER         gridID,
!                                     INTEGER         lpack)
      EXTERNAL        gridDefComplexPacking

      INTEGER         gridInqComplexPacking
!                                    (INTEGER         gridID)
      EXTERNAL        gridInqComplexPacking

!
!  ZAXIS routines
!
!                     zaxisName
!                                    (INTEGER         zaxistype,
!                                     CHARACTER*(*)   zaxisname)
      EXTERNAL        zaxisName

      INTEGER         zaxisCreate
!                                    (INTEGER         zaxistype,
!                                     INTEGER         size)
      EXTERNAL        zaxisCreate

!                     zaxisDestroy
!                                    (INTEGER         zaxisID)
      EXTERNAL        zaxisDestroy

      INTEGER         zaxisInqType
!                                    (INTEGER         zaxisID)
      EXTERNAL        zaxisInqType

      INTEGER         zaxisInqSize
!                                    (INTEGER         zaxisID)
      EXTERNAL        zaxisInqSize

      INTEGER         zaxisDuplicate
!                                    (INTEGER         zaxisID)
      EXTERNAL        zaxisDuplicate

!                     zaxisResize
!                                    (INTEGER         zaxisID,
!                                     INTEGER         size)
      EXTERNAL        zaxisResize

!                     zaxisPrint
!                                    (INTEGER         zaxisID,
!                                     INTEGER         index)
      EXTERNAL        zaxisPrint

!                     zaxisDefLevels
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION levels_vec)
      EXTERNAL        zaxisDefLevels

!                     zaxisInqLevels
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION levels_vec)
      EXTERNAL        zaxisInqLevels

!                     zaxisDefLevel
!                                    (INTEGER         zaxisID,
!                                     INTEGER         levelID,
!                                     DOUBLEPRECISION levels)
      EXTERNAL        zaxisDefLevel

      DOUBLEPRECISION zaxisInqLevel
!                                    (INTEGER         zaxisID,
!                                     INTEGER         levelID)
      EXTERNAL        zaxisInqLevel

!                     zaxisDefNlevRef
!                                    (INTEGER         gridID,
!                                     INTEGER         nhlev)
      EXTERNAL        zaxisDefNlevRef

      INTEGER         zaxisInqNlevRef
!                                    (INTEGER         gridID)
      EXTERNAL        zaxisInqNlevRef

!                     zaxisDefNumber
!                                    (INTEGER         gridID,
!                                     INTEGER         number)
      EXTERNAL        zaxisDefNumber

      INTEGER         zaxisInqNumber
!                                    (INTEGER         gridID)
      EXTERNAL        zaxisInqNumber

!                     zaxisDefUUID
!                                    (INTEGER         zaxisID,
!                                     INTEGER*1(16)   uuid)
      EXTERNAL        zaxisDefUUID

!                     zaxisInqUUID
!                                    (INTEGER         zaxisID,
!                                     INTEGER*1(16)   uuid)
      EXTERNAL        zaxisInqUUID

!                     zaxisDefName
!                                    (INTEGER         zaxisID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        zaxisDefName

!                     zaxisInqName
!                                    (INTEGER         zaxisID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        zaxisInqName

!                     zaxisDefLongname
!                                    (INTEGER         zaxisID,
!                                     CHARACTER*(*)   longname)
      EXTERNAL        zaxisDefLongname

!                     zaxisInqLongname
!                                    (INTEGER         zaxisID,
!                                     CHARACTER*(*)   longname)
      EXTERNAL        zaxisInqLongname

!                     zaxisDefUnits
!                                    (INTEGER         zaxisID,
!                                     CHARACTER*(*)   units)
      EXTERNAL        zaxisDefUnits

!                     zaxisInqUnits
!                                    (INTEGER         zaxisID,
!                                     CHARACTER*(*)   units)
      EXTERNAL        zaxisInqUnits

!                     zaxisInqStdname
!                                    (INTEGER         zaxisID,
!                                     CHARACTER*(*)   stdname)
      EXTERNAL        zaxisInqStdname

!                     zaxisDefPrec
!                                    (INTEGER         zaxisID,
!                                     INTEGER         prec)
      EXTERNAL        zaxisDefPrec

      INTEGER         zaxisInqPrec
!                                    (INTEGER         zaxisID)
      EXTERNAL        zaxisInqPrec

!                     zaxisDefPositive
!                                    (INTEGER         zaxisID,
!                                     INTEGER         positive)
      EXTERNAL        zaxisDefPositive

      INTEGER         zaxisInqPositive
!                                    (INTEGER         zaxisID)
      EXTERNAL        zaxisInqPositive

!                     zaxisDefLtype
!                                    (INTEGER         zaxisID,
!                                     INTEGER         ltype)
      EXTERNAL        zaxisDefLtype

      INTEGER         zaxisInqLtype
!                                    (INTEGER         zaxisID)
      EXTERNAL        zaxisInqLtype

!                     zaxisDefVct
!                                    (INTEGER         zaxisID,
!                                     INTEGER         size,
!                                     DOUBLEPRECISION vct_vec)
      EXTERNAL        zaxisDefVct

!                     zaxisInqVct
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION vct_vec)
      EXTERNAL        zaxisInqVct

      INTEGER         zaxisInqVctSize
!                                    (INTEGER         zaxisID)
      EXTERNAL        zaxisInqVctSize

!                     zaxisDefLbounds
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION lbounds_vec)
      EXTERNAL        zaxisDefLbounds

      INTEGER         zaxisInqLbounds
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION lbounds_vec)
      EXTERNAL        zaxisInqLbounds

      DOUBLEPRECISION zaxisInqLbound
!                                    (INTEGER         zaxisID,
!                                     INTEGER         index)
      EXTERNAL        zaxisInqLbound

!                     zaxisDefUbounds
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION ubounds_vec)
      EXTERNAL        zaxisDefUbounds

      INTEGER         zaxisInqUbounds
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION ubounds_vec)
      EXTERNAL        zaxisInqUbounds

      DOUBLEPRECISION zaxisInqUbound
!                                    (INTEGER         zaxisID,
!                                     INTEGER         index)
      EXTERNAL        zaxisInqUbound

!                     zaxisDefWeights
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION weights_vec)
      EXTERNAL        zaxisDefWeights

      INTEGER         zaxisInqWeights
!                                    (INTEGER         zaxisID,
!                                     DOUBLEPRECISION weights_vec)
      EXTERNAL        zaxisInqWeights

!                     zaxisChangeType
!                                    (INTEGER         zaxisID,
!                                     INTEGER         zaxistype)
      EXTERNAL        zaxisChangeType

!
!  TAXIS routines
!
      INTEGER         taxisCreate
!                                    (INTEGER         timetype)
      EXTERNAL        taxisCreate

!                     taxisDestroy
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisDestroy

      INTEGER         taxisDuplicate
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisDuplicate

!                     taxisCopyTimestep
!                                    (INTEGER         taxisIDdes,
!                                     INTEGER         taxisIDsrc)
      EXTERNAL        taxisCopyTimestep

!                     taxisDefType
!                                    (INTEGER         taxisID,
!                                     INTEGER         type)
      EXTERNAL        taxisDefType

!                     taxisDefVdate
!                                    (INTEGER         taxisID,
!                                     INTEGER         date)
      EXTERNAL        taxisDefVdate

!                     taxisDefVtime
!                                    (INTEGER         taxisID,
!                                     INTEGER         time)
      EXTERNAL        taxisDefVtime

      INTEGER         taxisInqVdate
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqVdate

      INTEGER         taxisInqVtime
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqVtime

!                     taxisDefRdate
!                                    (INTEGER         taxisID,
!                                     INTEGER         date)
      EXTERNAL        taxisDefRdate

!                     taxisDefRtime
!                                    (INTEGER         taxisID,
!                                     INTEGER         time)
      EXTERNAL        taxisDefRtime

      INTEGER         taxisInqRdate
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqRdate

      INTEGER         taxisInqRtime
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqRtime

!                     taxisDefFdate
!                                    (INTEGER         taxisID,
!                                     INTEGER         date)
      EXTERNAL        taxisDefFdate

!                     taxisDefFtime
!                                    (INTEGER         taxisID,
!                                     INTEGER         time)
      EXTERNAL        taxisDefFtime

      INTEGER         taxisInqFdate
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqFdate

      INTEGER         taxisInqFtime
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqFtime

      INTEGER         taxisHasBounds
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisHasBounds

!                     taxisDeleteBounds
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisDeleteBounds

!                     taxisDefVdateBounds
!                                    (INTEGER         taxisID,
!                                     INTEGER         vdate_lb,
!                                     INTEGER         vdate_ub)
      EXTERNAL        taxisDefVdateBounds

!                     taxisDefVtimeBounds
!                                    (INTEGER         taxisID,
!                                     INTEGER         vtime_lb,
!                                     INTEGER         vtime_ub)
      EXTERNAL        taxisDefVtimeBounds

!                     taxisInqVdateBounds
!                                    (INTEGER         taxisID,
!                                     INTEGER         vdate_lb,
!                                     INTEGER         vdate_ub)
      EXTERNAL        taxisInqVdateBounds

!                     taxisInqVtimeBounds
!                                    (INTEGER         taxisID,
!                                     INTEGER         vtime_lb,
!                                     INTEGER         vtime_ub)
      EXTERNAL        taxisInqVtimeBounds

!                     taxisDefCalendar
!                                    (INTEGER         taxisID,
!                                     INTEGER         calendar)
      EXTERNAL        taxisDefCalendar

      INTEGER         taxisInqCalendar
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqCalendar

!                     taxisDefTunit
!                                    (INTEGER         taxisID,
!                                     INTEGER         tunit)
      EXTERNAL        taxisDefTunit

      INTEGER         taxisInqTunit
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqTunit

!                     taxisDefForecastTunit
!                                    (INTEGER         taxisID,
!                                     INTEGER         tunit)
      EXTERNAL        taxisDefForecastTunit

      INTEGER         taxisInqForecastTunit
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqForecastTunit

!                     taxisDefForecastPeriod
!                                    (INTEGER         taxisID,
!                                     DOUBLEPRECISION fc_period)
      EXTERNAL        taxisDefForecastPeriod

      DOUBLEPRECISION taxisInqForecastPeriod
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqForecastPeriod

!                     taxisDefNumavg
!                                    (INTEGER         taxisID,
!                                     INTEGER         numavg)
      EXTERNAL        taxisDefNumavg

      INTEGER         taxisInqType
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqType

      INTEGER         taxisInqNumavg
!                                    (INTEGER         taxisID)
      EXTERNAL        taxisInqNumavg

      CHARACTER(80)   tunitNamePtr
!                                    (INTEGER         tunitID)
      EXTERNAL        tunitNamePtr

!
!  Institut routines
!
      INTEGER         institutDef
!                                    (INTEGER         center,
!                                     INTEGER         subcenter,
!                                     CHARACTER*(*)   name,
!                                     CHARACTER*(*)   longname)
      EXTERNAL        institutDef

      INTEGER         institutInq
!                                    (INTEGER         center,
!                                     INTEGER         subcenter,
!                                     CHARACTER*(*)   name,
!                                     CHARACTER*(*)   longname)
      EXTERNAL        institutInq

      INTEGER         institutInqNumber
      EXTERNAL        institutInqNumber

      INTEGER         institutInqCenter
!                                    (INTEGER         instID)
      EXTERNAL        institutInqCenter

      INTEGER         institutInqSubcenter
!                                    (INTEGER         instID)
      EXTERNAL        institutInqSubcenter

      CHARACTER(80)   institutInqNamePtr
!                                    (INTEGER         instID)
      EXTERNAL        institutInqNamePtr

      CHARACTER(80)   institutInqLongnamePtr
!                                    (INTEGER         instID)
      EXTERNAL        institutInqLongnamePtr

!
!  Model routines
!
      INTEGER         modelDef
!                                    (INTEGER         instID,
!                                     INTEGER         modelgribID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        modelDef

      INTEGER         modelInq
!                                    (INTEGER         instID,
!                                     INTEGER         modelgribID,
!                                     CHARACTER*(*)   name)
      EXTERNAL        modelInq

      INTEGER         modelInqInstitut
!                                    (INTEGER         modelID)
      EXTERNAL        modelInqInstitut

      INTEGER         modelInqGribID
!                                    (INTEGER         modelID)
      EXTERNAL        modelInqGribID

      CHARACTER(80)   modelInqNamePtr
!                                    (INTEGER         modelID)
      EXTERNAL        modelInqNamePtr

!
!  Table routines
!
!                     tableWriteC
!                                    (CHARACTER*(*)   filename,
!                                     INTEGER         tableID)
      EXTERNAL        tableWriteC

!                     tableWrite
!                                    (CHARACTER*(*)   filename,
!                                     INTEGER         tableID)
      EXTERNAL        tableWrite

      INTEGER         tableRead
!                                    (CHARACTER*(*)   tablefile)
      EXTERNAL        tableRead

      INTEGER         tableDef
!                                    (INTEGER         modelID,
!                                     INTEGER         tablenum,
!                                     CHARACTER*(*)   tablename)
      EXTERNAL        tableDef

      CHARACTER(80)   tableInqNamePtr
!                                    (INTEGER         tableID)
      EXTERNAL        tableInqNamePtr

!                     tableDefEntry
!                                    (INTEGER         tableID,
!                                     INTEGER         code,
!                                     CHARACTER*(*)   name,
!                                     CHARACTER*(*)   longname,
!                                     CHARACTER*(*)   units)
      EXTERNAL        tableDefEntry

      INTEGER         tableInq
!                                    (INTEGER         modelID,
!                                     INTEGER         tablenum,
!                                     CHARACTER*(*)   tablename)
      EXTERNAL        tableInq

      INTEGER         tableInqNumber
      EXTERNAL        tableInqNumber

      INTEGER         tableInqNum
!                                    (INTEGER         tableID)
      EXTERNAL        tableInqNum

      INTEGER         tableInqModel
!                                    (INTEGER         tableID)
      EXTERNAL        tableInqModel

!                     tableInqPar
!                                    (INTEGER         tableID,
!                                     INTEGER         code,
!                                     CHARACTER*(*)   name,
!                                     CHARACTER*(*)   longname,
!                                     CHARACTER*(*)   units)
      EXTERNAL        tableInqPar

      INTEGER         tableInqParCode
!                                    (INTEGER         tableID,
!                                     CHARACTER*(*)   name,
!                                     INTEGER         code)
      EXTERNAL        tableInqParCode

      INTEGER         tableInqParName
!                                    (INTEGER         tableID,
!                                     INTEGER         code,
!                                     CHARACTER*(*)   name)
      EXTERNAL        tableInqParName

      INTEGER         tableInqParLongname
!                                    (INTEGER         tableID,
!                                     INTEGER         code,
!                                     CHARACTER*(*)   longname)
      EXTERNAL        tableInqParLongname

      INTEGER         tableInqParUnits
!                                    (INTEGER         tableID,
!                                     INTEGER         code,
!                                     CHARACTER*(*)   units)
      EXTERNAL        tableInqParUnits

      CHARACTER(80)   tableInqParNamePtr
!                                    (INTEGER         tableID,
!                                     INTEGER         parID)
      EXTERNAL        tableInqParNamePtr

      CHARACTER(80)   tableInqParLongnamePtr
!                                    (INTEGER         tableID,
!                                     INTEGER         parID)
      EXTERNAL        tableInqParLongnamePtr

      CHARACTER(80)   tableInqParUnitsPtr
!                                    (INTEGER         tableID,
!                                     INTEGER         parID)
      EXTERNAL        tableInqParUnitsPtr

!
!  History routines
!
!                     streamDefHistory
!                                    (INTEGER         streamID,
!                                     INTEGER         size,
!                                     CHARACTER*(*)   history)
      EXTERNAL        streamDefHistory

      INTEGER         streamInqHistorySize
!                                    (INTEGER         streamID)
      EXTERNAL        streamInqHistorySize

!                     streamInqHistoryString
!                                    (INTEGER         streamID,
!                                     CHARACTER*(*)   history)
      EXTERNAL        streamInqHistoryString

!
!  Subtype routines
!
      INTEGER         subtypeCreate
!                                    (INTEGER         subtype)
      EXTERNAL        subtypeCreate

!
!  Gives a textual summary of the variable subtype
!
!                     subtypePrint
!                                    (INTEGER         subtypeID)
      EXTERNAL        subtypePrint

!
!  Compares two subtype data structures.
!
      INTEGER         subtypeCompare
!                                    (INTEGER         subtypeID1,
!                                     INTEGER         subtypeID2)
      EXTERNAL        subtypeCompare

      INTEGER         subtypeInqSize
!                                    (INTEGER         subtypeID)
      EXTERNAL        subtypeInqSize

      INTEGER         subtypeInqActiveIndex
!                                    (INTEGER         subtypeID)
      EXTERNAL        subtypeInqActiveIndex

!                     subtypeDefActiveIndex
!                                    (INTEGER         subtypeID,
!                                     INTEGER         index)
      EXTERNAL        subtypeDefActiveIndex

!
!  Generate a "query object" out of a key-value pair.
!
!
!  Generate an AND-combined "query object" out of two previous
!  query objects.
!
      INTEGER         subtypeInqTile
!                                    (INTEGER         subtypeID,
!                                     INTEGER         tileindex,
!                                     INTEGER         attribute)
      EXTERNAL        subtypeInqTile

      INTEGER         vlistInqVarSubtype
!                                    (INTEGER         vlistID,
!                                     INTEGER         varID)
      EXTERNAL        vlistInqVarSubtype

!                     gribapiLibraryVersion
!                                    (INTEGER         major_version,
!                                     INTEGER         minor_version,
!                                     INTEGER         revision_version)
      EXTERNAL        gribapiLibraryVersion

