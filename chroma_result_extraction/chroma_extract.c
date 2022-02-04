//-------------------------------------------------------------------------------------------//
//-- Kim Somfleth <kim.somfleth@internode.on.net>                                          --//
//--                        LIME FILE EXTRACTER FOR FEYNMAN HELLMANN                       --//
//-- 01/02/16 CREATED STUB                                                                 --//
//-------------------------------------------------------------------------------------------//
//. build_lime.sh && (cd ../../build/harmless && ./hl_lime 0 21 hadron_qcdsf* ); cd ../../src/harmless

#include "../hl.h"
#include "chroma_lime_reader.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "chroma_momentum_numbers.h"

#define NUM_CORRELATOR_SETS 65536 //WE DO THIS BECAUSE WE TRY A STATIC HASH MAP RATHER THAN A LINKED LIST ONE FOR SIMPLICITY 
// FROM /home/566/kys566/opt/src/jz_chroma/lib/meas/hadron/barhqlq_qcdsf_w.cc
#if 0 
#define NUCLEON_UP    0  //-- ncln1        --//
#define NUCLEON_DOWN 21  //-- ncln2        --//
#define DELTA_3HALF   2  //-- delta_1      --//
#define DELTA_M3HALF 25  //-- delta_m3half --//
#define DELTA_HALF   23  //-- delta_half   --//
#define DELTA_MHALF  24  //-- delta_mhalf  --//
#define PION_G5G5    255 //-- pion_g5_g5   --//
#endif
c8 *baryon_identity_string[] = {
  "nucleon"    //0
, ""           //1
, "delta3half" //2
, ""           //3
, ""           //4
, ""           //5
, ""           //6
, ""           //7
, ""           //8
, "nucleon"    //9
, ""           //10
, ""           //11
, ""           //12
, ""           //13
, ""           //14
, ""           //15
, ""           //16
, "nucleon"    //17
, "nucleon"    //18
, "nucleon"    //19
, ""           //20
, "nucleon"    //21
, ""           //22
, "deltahalf"  //23
, "deltahalf"  //24
, "delta3half" //25
};
c8 *spin_identity_string[] = {
  "up"     //0
, ""       //1
, "up"     //2
, ""       //3
, ""       //4
, ""       //5
, ""       //6
, ""       //7
, ""       //8
, "unpol"  //9
, ""       //10
, ""       //11
, ""       //12
, ""       //13
, ""       //14
, ""       //15
, ""       //16
, "polx"   //17
, "poly"   //18
, "polz"   //19
, ""       //20
, "down"   //21
, ""       //22
, "up"     //23
, "down"   //24
, "down"   //25
};

c8 *real_operator_identity_string[] = {
  "I"    //0
, "g1"   //1
, "g2"   //2
, "g1g2" //3
, "g3"   //4
, "g1g3" //5
, "g2g3" //6
, "g5g4" //7
, "g4"   //8
, "g1g4" //9
, "g2g4" //10
, "g3g5" //11
, "g3g4" //12
, "g6g2" //13
, "g1g5" //14
, "g5"   //15
};
c8 *imag_operator_identity_string[] = {
  "iI"    //0
, "ig1"   //1
, "ig2"   //2
, "ig1g2" //3
, "ig3"   //4
, "ig1g3" //5
, "ig2g3" //6
, "ig5g4" //7
, "ig4"   //8
, "ig1g4" //9
, "ig2g4" //10
, "ig3g5" //11
, "ig3g4" //12
, "ig6g2" //13
, "ig1g5" //14
, "ig5g5" //15
  "cv1"   //-- CONSERVED CURRENTS HAVE JUST ADDED A BIT TO THE OPCODE --//
, "cv2"
, "cv3"
, "cv4"
};

u64 write( i32 dest, void *buffer, u64 buffer_size ); //-- INCLUDES UNISTR --//

typedef enum {
  hlSHELL_SMEARING = 1
, hlPOINT_SMEARING
} hlSmearing;

c8 *hlSmearingStrs[] = {
  "undefined_smear"
, "sh"
, "pt"
};

typedef struct Gauge_Chunk {
  struct Gauge_Chunk *next;
         c8          *ident_string;
         u64          source[ 4 ];
         u64          r64_count;
         r64         *data;
         r64         *data_trev;
         r64         *forward_test_correlator[3];
} Gauge_Chunk;

typedef enum {
  hlBARYON
, hlMESON
} hlHadron;

c8 *hlHadronStrs[] = {
  "Baryon"
, "Meson"
};

typedef struct Correlator {
  struct Correlator *next;
         c8         *ident_string;
         u64         source[4];
         u64         r64_count;
         r64        *data;
         r64        *data_trev;
         r64        *forward_test_correlator[3];
} Correlator;

typedef struct Feynman_Hellmann_Params {
  r64 lambda[2];
  i64 chroma_opcode;
  b32 conserved_current;
  i64 lattice_momentum_transfer[3];
} Feynman_Hellmann_Params;

typedef struct newPropagator {
  r64 prop_kappa;
  hlSmearing sync_smear;
  hlSmearing src_smear;
  u64        lambda_count;
  Feynman_Hellmann_Params feynman_hellmann_params[16];
} newPropagator;

typedef struct Correlator_Set {
  u64         hash;
  hlHadron    hadron_type;
  u64         dim[4];
  u64         hadron_number;
  r64         sea_kappa1;
  r64         sea_kappa2;
  i64         lattice_momentum[3];
  u64         prop_count;
  newPropagator props[3];
  u64         bad_correlators;
  u64         correlator_count;
  Correlator *first_correlator;
  Correlator *last_correlator;
} Correlator_Set;

static void PrintCorrelatorSet( Correlator_Set cor_set ){
  printf( "-------------------------\n");
  printf( "CorrelatorSet: %s (%lu)\n", hlHadronStrs[cor_set.hadron_type], cor_set.hash );
  printf( "Dim: (%lu,%lu,%lu,%lu)\n", cor_set.dim[0], cor_set.dim[1], cor_set.dim[2], cor_set.dim[3] );
  printf( "Hadron Numbers %lu\n", cor_set.hadron_number );
  printf( "Sea Kappas (%f,%f)\n", cor_set.sea_kappa1, cor_set.sea_kappa2 );
  printf( "Sync Mom (%li,%li,%li)\n", cor_set.lattice_momentum[0], cor_set.lattice_momentum[1], cor_set.lattice_momentum[2] );
  printf( "Extraction of %lu good cors (%lu bad)\n", cor_set.correlator_count, cor_set.bad_correlators );
  for( u64 prop_index = 0; prop_index < cor_set.prop_count; ++prop_index ){
    newPropagator *prop = cor_set.props + prop_index;
    printf( "\tProp Kappa %f\n", prop->prop_kappa );
    printf( "\tProp smearing sync: %s, src: %s\n", hlSmearingStrs[prop->sync_smear], hlSmearingStrs[prop->src_smear] );
    for( u64 lambda_index = 0; lambda_index < prop->lambda_count; ++lambda_index ){
      Feynman_Hellmann_Params *fhp = prop->feynman_hellmann_params + lambda_index;
      printf( "\t\t (%f,%f), %s, q:(%li,%li,%li)\n", fhp->lambda[0], fhp->lambda[1], real_operator_identity_string[fhp->chroma_opcode], fhp->lattice_momentum_transfer[0], fhp->lattice_momentum_transfer[1], fhp->lattice_momentum_transfer[2] );
    }
  }
}

static hlSmearing GetSmearingType( c8 *str ){
  hlSmearing result = 0;
  if( hlStrStartsWithC( "POINT_", str ) ){
    result = hlPOINT_SMEARING;
  }else if( hlStrStartsWithC( "SHELL_", str ) ){
    result = hlSHELL_SMEARING;
  }
  hlASSERT( result );
  //hlASSERT( result == hlSHELL_SMEARING || result == hlPOINT_SMEARING );
  return( result );
}

static inline u64 HashU64(u64 x) {
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
  x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
  x = x ^ (x >> 31);
  return x;
}

static u64 HashCorrelatorSet( Correlator_Set *set ){
  u64 hash = 0;
  hash += HashU64( 1922*set->hadron_type + 12823*set->hadron_number + (u64)(1048*(set->sea_kappa1+set->sea_kappa2)) + 12283928312*set->dim[0] + 123878*set->dim[1] + 12938998*set->dim[2] + 982139128*set->dim[3] + 12391289*set->lattice_momentum[0] + 92139812*set->lattice_momentum[1] + 9812319*set->lattice_momentum[2] );
  for( u64 prop_index = 0; prop_index < set->prop_count; ++prop_index ){
    newPropagator *prop = set->props + prop_index;
    hash += HashU64( (u64)(2049*prop->prop_kappa) + 123909*prop->sync_smear + 84520*prop->src_smear + 8*prop_index );
    for( u64 lambda_index = 0; lambda_index < prop->lambda_count; ++lambda_index ){
      Feynman_Hellmann_Params *params = prop->feynman_hellmann_params + lambda_index;
      hash += HashU64( ( lambda_index + 8 ) * ( lambda_index + prop_index + 4 ) * ( (u64)(1028*params->lambda[0] + 10923*params->lambda[1]) + 1239812*params->chroma_opcode + 81238923*params->lattice_momentum_transfer[0] + 8402083*params->lattice_momentum_transfer[1] + 12398129*params->lattice_momentum_transfer[2] ) );
    }
  }
  return( hash );
}

static b32 CorrelatorSetCmp( Correlator_Set *set1, Correlator_Set *set2 ){
  b32 same_hadron_type        = set1->hadron_type   == set2->hadron_type;
  b32 same_hadron_number      = set1->hadron_number == set2->hadron_number;
  b32 same_sea_kappa1         = set1->sea_kappa1    == set2->sea_kappa1;
  b32 same_sea_kappa2         = set1->sea_kappa2    == set2->sea_kappa2;
  b32 same_prop_count         = set1->prop_count    == set2->prop_count;
  b32 same_dim_0              = set1->dim[0]        == set2->dim[0];
  b32 same_dim_1              = set1->dim[1]        == set2->dim[1];
  b32 same_dim_2              = set1->dim[2]        == set2->dim[2];
  b32 same_dim_3              = set1->dim[3]        == set2->dim[3];
  b32 same_lattice_momentum_0 = set1->lattice_momentum[0] == set2->lattice_momentum[0];
  b32 same_lattice_momentum_1 = set1->lattice_momentum[1] == set2->lattice_momentum[1];
  b32 same_lattice_momentum_2 = set1->lattice_momentum[2] == set2->lattice_momentum[2];
  b32 same_prop               = set1->prop_count          == set2->prop_count;
  for( u64 prop_index = 0; same_prop && prop_index < set1->prop_count; ++prop_index ){
    b32 same_prop_kappa  = set1->props[prop_index].prop_kappa   == set2->props[prop_index].prop_kappa;
    b32 same_sync_smear  = set1->props[prop_index].sync_smear   == set2->props[prop_index].sync_smear; 
    b32 same_src_smear   = set1->props[prop_index].src_smear    == set2->props[prop_index].src_smear; 
    b32 same_fhp         = set1->props[prop_index].lambda_count == set2->props[prop_index].lambda_count;
    for( u64 lambda_index = 0; same_fhp && lambda_index < set1->props[prop_index].lambda_count; ++lambda_index ){
      b32 same_lambda0 = set1->props[prop_index].feynman_hellmann_params[lambda_index].lambda[0]
                      ==  set2->props[prop_index].feynman_hellmann_params[lambda_index].lambda[0];
      b32 same_lambda1 = set1->props[prop_index].feynman_hellmann_params[lambda_index].lambda[1]
                      ==  set2->props[prop_index].feynman_hellmann_params[lambda_index].lambda[1];
      b32 same_op      = set1->props[prop_index].feynman_hellmann_params[lambda_index].chroma_opcode
                      ==  set2->props[prop_index].feynman_hellmann_params[lambda_index].chroma_opcode;
      b32 same_q_0     = set1->props[prop_index].feynman_hellmann_params[lambda_index].lattice_momentum_transfer[0]
                      ==  set2->props[prop_index].feynman_hellmann_params[lambda_index].lattice_momentum_transfer[0];
      b32 same_q_1     = set1->props[prop_index].feynman_hellmann_params[lambda_index].lattice_momentum_transfer[1]
                      ==  set2->props[prop_index].feynman_hellmann_params[lambda_index].lattice_momentum_transfer[1];
      b32 same_q_2     = set1->props[prop_index].feynman_hellmann_params[lambda_index].lattice_momentum_transfer[2]
                      ==  set2->props[prop_index].feynman_hellmann_params[lambda_index].lattice_momentum_transfer[2];
      same_fhp &= same_lambda0 && same_lambda1 && same_op && same_q_0 && same_q_1 && same_q_2;
    }
    same_prop &= same_prop_kappa && same_sync_smear && same_src_smear && same_fhp;
  }
  b32 result = same_hadron_type && same_hadron_number && same_sea_kappa1 && same_sea_kappa2 && same_prop_count && same_dim_0 && same_dim_1 && same_dim_2 && same_dim_3 && same_lattice_momentum_0 && same_lattice_momentum_1 && same_lattice_momentum_2 && same_prop;

  hlASSERT( HashCorrelatorSet( set1 ) == HashCorrelatorSet( set2 ) );
  return( result );
}

static Correlator_Set *GetOrCopyCorrelatorSet( hlHash_Map *hash_map, Correlator_Set *correlator_sets, Correlator_Set params ){
  u64 hash = HashCorrelatorSet( &params );
  //printf( "[DEBUG] FINDING HASH %lu\n", hash );
  Correlator_Set *result = hlHashMapGet( hash_map, hash, Correlator_Set );
  if( !result ){
    //hlDEBUG_STR( "WASN'T THERE GENERATING\n" );
    while( correlator_sets->hash ){
      ++correlator_sets;
    }
    hlMemCopy( &params, correlator_sets, sizeof(params) );
    result = correlator_sets;
    result->hash = hash;
    hlHashMapInsert( hash_map, result->hash, result );
  }else{
    //hlDEBUG_STR( "IS ALREADY THERE\n" );
    if( !CorrelatorSetCmp( result, &params ) ){
      hlDEBUG_STR( "TESTING WHETHER THE HASHER ISN'T CONFLICTING\n" );
      PrintCorrelatorSet( *result );
      params.hash = hash;
      PrintCorrelatorSet( params );
    }
    hlASSERT( CorrelatorSetCmp( result, &params ) );
  }
  return( result );
}

static hlXML_Result SpecialStupidChromaForwardTestPropQuery( hlMem_Pool *tmp_pool, hlMem_Pool *pool, hlXML *root, c8 *query ){
  hlXML_Result result = hlXMLQuery( tmp_pool, pool, root, query );
  if( result.count == 1 ){
    hlXML *xml = result.values[ 0 ];
    if( xml->type == hlXML_NUM_ARR ){
      //-- BECAUSE FLOATS ARE ACTUALLY INTS IF YOU SQUINT HARD ENOUGH OF COURSE THIS ONLY OCCURS IN EVERY 1000 TRAJECTORIES --//
      r64 *the_actual_values = hlPUSHARRAY( pool, xml->num_arr_val.count, r64 );
      for( u64 number_index = 0; number_index < xml->num_arr_val.count; ++number_index ){
        the_actual_values[ number_index ] = (r64)xml->num_arr_val.values[ number_index ];
      }
      xml->type               = hlXML_FLT_ARR;
      xml->flt_arr_val.count  = xml->num_arr_val.count;
      xml->flt_arr_val.values = the_actual_values;
    }
  }
  return( result );
}

i32 main( i32 argc, c8 **argv ){
  hlMem_Pool mem_pool;
  hlMem_Pool tmp_pool;
  hlInitializeProgramMemory( &mem_pool, GIGABYTES( 32 ), &tmp_pool, GIGABYTES( 6 ) );

  u64 hadron_numbers_count = argc - 2;
  hlASSERT( hadron_numbers_count > 0 );
  u64 *hadron_numbers = hlPUSHARRAY( &mem_pool, hadron_numbers_count, u64 );
  for( u64 arg_index = 1; arg_index < hadron_numbers_count + 1; ++arg_index ){
    hadron_numbers[ arg_index - 1 ] = atol( argv[ arg_index ] );
  }

  c8 *lime_listing_file_name = argv[ argc-1 ];
  hlFile lime_listing = hlFileRead( &mem_pool, lime_listing_file_name );
  u64 lime_files = 0;
  c8 **lime_filenames = hlStrSplitC( &mem_pool, '\n', &lime_files, (c8 *)lime_listing.contents );

  //Correlator_Set correlator_sets[NUM_CORRELATOR_SETS] = {};
  Correlator_Set *correlator_sets = hlPUSHARRAY( &mem_pool, NUM_CORRELATOR_SETS, Correlator_Set );
  hlHash_Map correlator_hashmap = hlHashMapGenerate( &mem_pool, 3*NUM_CORRELATOR_SETS );


  printf( "PARSING: %lu FILES\n", lime_files );
  for( u64 lime_index = 0; lime_index < lime_files; ++lime_index ){
    c8   *path = lime_filenames[ lime_index ];
    FILE *file = fopen( path, "r" );
    hlASSERT( file );
    
    fseek( file, 0, SEEK_END );
    u64 file_size = ftell( file );
    fseek( file, 0, SEEK_SET );

    hlResetMemoryPool( &tmp_pool );

    u8 *file_contents = hlPUSHARRAY( &tmp_pool, file_size, u8 );
    fread( file_contents, file_size, 1, file );
    fclose( file );

    u8 *current_position = file_contents;
    
    Lime_Header *lime_header   = ReadLimeHeader( &current_position );
    hlASSERT( hlStrStartsWithC( "qcdsfDir", lime_header->type ) );

    hlXML dir_xml     = hlParseXML( &tmp_pool, &mem_pool, (c8 *)lime_header->data, lime_header->data_size );
    hlXML *sink_pairs = hlXMLQueryArr(    &tmp_pool, &tmp_pool, &dir_xml, "Input/NamedObject/sink_pairs"     );
    i64 *dim          = hlXMLQueryNumVec( &tmp_pool, &tmp_pool, &dir_xml, "ProgramInfo/Setgeom/latt_size", 4 );
    i64  mom2_max     = hlXMLQueryNum(    &tmp_pool, &tmp_pool, &dir_xml, "Input/Param/mom2_max"             );
    //hlXML_Result configuration_xml   = hlXMLQuery( &tmp_pool, &tmp_pool, &dir_xml, "Config_info/ILDG/LFN"          );
    //hlASSERT( configuration_xml.count == 1                       );
    //hlASSERT( configuration_xml.values[ 0 ]->type == hlXML_STR   );

    //hlDEBUG_STR( "Getting config name\n" );
    c8 *ident_str = GetConfigurationName( &tmp_pool, &mem_pool, &dir_xml, path );
    //if( !configuration_xml.count ){
    //  //hlDEBUG_STR( "Alternate Lime Format\n" );
    //  configuration_xml = hlXMLQuery( &tmp_pool, &tmp_pool, &dir_xml, "Input/lime_file" );
    //  u64 tail_index = 1 + hlStrFindLastOccurenceC( configuration_xml.values[ 0 ]->str, '/' );
    //  tail_index += 1 + hlStrFindFirstOccurenceC( configuration_xml.values[ 0 ]->str + tail_index, '_' );
    //  u64 last_index = hlStrFindNthToLastOccurenceC( configuration_xml.values[ 0 ]->str + tail_index, '_', 6 );
    //  ident_str = hlStrCopy( &mem_pool, configuration_xml.values[ 0 ]->str + tail_index, last_index );
    //}else{
    //  u64 tail_index = 1 + hlStrFindLastOccurenceC( configuration_xml.values[ 0 ]->str, '/' );
    //  ident_str  = hlStrCopyC( &mem_pool, configuration_xml.values[ 0 ]->str + tail_index );
    //}
    //hlXMLPrintChildren( &dir_xml );
    printf( "PARSING: %s <- %s\n", ident_str, path );
    //PrintGaugeSet( gauge_set );

    u64 nt = dim[3];
    b32 first_bad_prop = 1;
    for( u64 prop_pair_index = 0; prop_pair_index < sink_pairs->obj_val.count; ++prop_pair_index ){
      Correlator_Set correlator_params = {};
      for( u64 dim_index = 0; dim_index < hlARRAYCOUNT( correlator_params.dim ); ++dim_index ){
        correlator_params.dim[dim_index] = dim[dim_index];
      }
      Correlator correlator = {};
      //printf( "PARSING META HEADER\n" );
      Lime_Header *meta_xml_header  = ReadLimeHeader( &current_position );
      hlASSERT( hlStrStartsWithC( "meta-xml", meta_xml_header->type ) );

      //printf( "PARSING HADRON BIN HEADER\n" );
      Lime_Header *hadron_bin_header = ReadLimeHeader( &current_position );
      b32 is_baryon_bin = hlStrStartsWithC( "baryons-bin", hadron_bin_header->type );
      b32 is_meson_bin  = hlStrStartsWithC( "mesons-bin" , hadron_bin_header->type );
      b32 is_baryon     = is_baryon_bin;
      correlator_params.hadron_type = is_baryon ? hlBARYON : hlMESON;
      hlASSERT( is_baryon_bin || is_meson_bin );
      hlASSERT( is_baryon_bin != is_meson_bin ); //-- RATHER THAN USE XOR TO WORK OUT PROBLEM IMMEDIATELY --//

      Lime_Header *hadron_trev_bin_header;
      if( is_baryon_bin ){ 
        //printf( "PARSING TREV HADRON BIN HEADER\n" );
        hadron_trev_bin_header = ReadLimeHeader( &current_position );
        b32 is_baryon_rev_bin = hlStrStartsWithC( "baryons-trev-bin", hadron_trev_bin_header->type );
        hlASSERT( is_baryon_rev_bin );
      }


#define PROP1_PATH              "Forward_prop_headers/First_forward_prop" 
#define PROP2_PATH              "Forward_prop_headers/Second_forward_prop"
#define PROP3_PATH              "Forward_prop_headers/Third_forward_prop" 
#define FORWARD_PROP_COR_PATH_1 "Forward_prop_correlator/forward_prop_corr_1"
#define FORWARD_PROP_COR_PATH_2 "Forward_prop_correlator/forward_prop_corr_2"
#define FORWARD_PROP_COR_PATH_3 "Forward_prop_correlator/forward_prop_corr_3"
      hlXML       meta_xml   = hlParseXML( &tmp_pool, &mem_pool, (c8 *)meta_xml_header->data, meta_xml_header->data_size );
      //printf( "%s\n", meta_xml_header->data );
      r64 mass_1 = hlXMLQueryFlt( &tmp_pool, &mem_pool, &meta_xml, "Mass_1" );
      r64 mass_2 = hlXMLQueryFlt( &tmp_pool, &mem_pool, &meta_xml, "Mass_2" );
      hlXML_Result prop_query_xml[ 3 ];
                   prop_query_xml[ 0 ] = hlXMLQuery( &tmp_pool, &mem_pool, &meta_xml, PROP1_PATH );
                   prop_query_xml[ 1 ] = hlXMLQuery( &tmp_pool, &mem_pool, &meta_xml, PROP2_PATH );
                   prop_query_xml[ 2 ] = hlXMLQuery( &tmp_pool, &mem_pool, &meta_xml, PROP3_PATH );
      hlXML_Result fwd_prop_cor_xml[ 3 ];
		   fwd_prop_cor_xml[ 0 ] = SpecialStupidChromaForwardTestPropQuery( &tmp_pool, &mem_pool, &meta_xml, FORWARD_PROP_COR_PATH_1 );
                   fwd_prop_cor_xml[ 1 ] = SpecialStupidChromaForwardTestPropQuery( &tmp_pool, &mem_pool, &meta_xml, FORWARD_PROP_COR_PATH_2 );
                   fwd_prop_cor_xml[ 2 ] = SpecialStupidChromaForwardTestPropQuery( &tmp_pool, &mem_pool, &meta_xml, FORWARD_PROP_COR_PATH_3 );
      u64 prop_count = 2;
      hlASSERT( prop_query_xml[ 0 ].count   == 1 );
      hlASSERT( prop_query_xml[ 1 ].count   == 1 );
      hlASSERT( fwd_prop_cor_xml[ 0 ].count == 1 );
      hlASSERT( fwd_prop_cor_xml[ 1 ].count == 1 );
      //hlXML *fwd_prop_cor = hlXMLQueryObject(  &tmp_pool, &mem_pool, &meta_xml, "Forward_prop_correlator" );
      //hlXMLPrintChildren( fwd_prop_cor );
      hlASSERT( fwd_prop_cor_xml[ 0 ].values[ 0 ]->type == hlXML_FLT_ARR );
      hlASSERT( fwd_prop_cor_xml[ 1 ].values[ 0 ]->type == hlXML_FLT_ARR );
      hlASSERT( fwd_prop_cor_xml[ 0 ].values[ 0 ]->flt_arr_val.count == nt );
      hlASSERT( fwd_prop_cor_xml[ 1 ].values[ 0 ]->flt_arr_val.count == nt );
      if( prop_query_xml[ 2 ].count ){
        ++prop_count;
        hlASSERT( prop_query_xml[ 2 ].count   == 1 );
        hlASSERT( fwd_prop_cor_xml[ 2 ].count == 1 );
        hlASSERT( fwd_prop_cor_xml[ 2 ].values[ 0 ]->type == hlXML_FLT_ARR );
        hlASSERT( fwd_prop_cor_xml[ 2 ].values[ 0 ]->flt_arr_val.count == nt );
        hlASSERT( !"Need to implement smearing query" );
      }
      r64 kappa1 = 1.0 / ( 2.0 * mass_1 + 8.0 );
      r64 kappa2 = 1.0 / ( 2.0 * mass_2 + 8.0 );
      correlator_params.sea_kappa1 = kappa1;
      correlator_params.sea_kappa2 = kappa2;
      correlator_params.prop_count = prop_count;

#define PROP_SRC_PATH                      "PropSource/Source/t_srce"
#define PROP_SRC_SMEAR_PATH                "PropSource/Source/SourceType"
#define PROP_SYNC_SMEAR_PATH               "PropSink/Sink/SinkType"
#define PROP_KAPPA_PATH                    "ForwardProp/FermionAction/Kappa"
      for( u64 prop_index = 0; prop_index < prop_count; ++prop_index ){
        hlXML *prop_xml = prop_query_xml[ prop_index ].values[ 0 ];
        hlASSERT( prop_xml );
  
        i64 *prop_source     = hlXMLQueryNumVec( &tmp_pool, &mem_pool, prop_xml, PROP_SRC_PATH                 , hlARRAYCOUNT(correlator.source) );
        c8  *prop_sync_smear = hlXMLQueryStr(    &tmp_pool, &mem_pool, prop_xml, PROP_SYNC_SMEAR_PATH           );
        c8  *prop_src_smear  = hlXMLQueryStr(    &tmp_pool, &mem_pool, prop_xml, PROP_SRC_SMEAR_PATH            );
        r64  prop_kappa      = hlXMLQueryFlt(    &tmp_pool, &mem_pool, prop_xml, PROP_KAPPA_PATH                );
        hlSmearing sync_smearing  = GetSmearingType( prop_sync_smear );
        hlSmearing src_smearing   = GetSmearingType( prop_src_smear  );
        for( u64 source_index = 0; source_index < hlARRAYCOUNT(correlator.source); ++source_index ){
          correlator.source[source_index] = prop_source[source_index];
        }
        correlator_params.props[prop_index].prop_kappa   = prop_kappa;
        correlator_params.props[prop_index].sync_smear   = sync_smearing;
        correlator_params.props[prop_index].src_smear    = src_smearing;
        correlator_params.props[prop_index].lambda_count = 0;

#define PROP_FEYNHELL_LAMBDA_REAL_PATH            "ForwardProp/FermionAction/FeynHellParam/LambdaReal"
#define PROP_FEYNHELL_LAMBDA_IMAG_PATH            "ForwardProp/FermionAction/FeynHellParam/LambdaImag"
#define PROP_FEYNHELL_OPERATOR_PATH               "ForwardProp/FermionAction/FeynHellParam/Operator"
#define PROP_FEYNHELL_MOMENTUM_PATH               "ForwardProp/FermionAction/FeynHellParam/Momentum"
#define PROP_FEYNHELL_LAMBDA_REAL_PATH_OLD        "ForwardProp/FermionAction/LambdaReal"
#define PROP_FEYNHELL_LAMBDA_IMAG_PATH_OLD        "ForwardProp/FermionAction/LambdaImag"
#define PROP_FEYNHELL_OPERATOR_PATH_OLD           "ForwardProp/FermionAction/Operator"
#define PROP_FEYNHELL_MOMENTUM_PATH_OLD           "ForwardProp/FermionAction/Momentum"
#define PROP_FEYNHELL_SOURCE_MOD_LAMBDA_REAL_PATH "ForwardProp/FermionAction/FermState/FeynHellParam/LambdaReal"
#define PROP_FEYNHELL_SOURCE_MOD_LAMBDA_IMAG_PATH "ForwardProp/FermionAction/FermState/FeynHellParam/LambdaImag"
#define PROP_FEYNHELL_SOURCE_MOD_OPERATOR_PATH    "ForwardProp/FermionAction/FermState/FeynHellParam/Operator"
#define PROP_FEYNHELL_SOURCE_MOD_MOMENTUM_PATH    "ForwardProp/FermionAction/FermState/FeynHellParam/Momentum"
        typedef enum hl_feynhell_type {
          hlCHROMA_OLD_SINGLE_INVERSION_MOD
        , hlCHROMA_INVERSION_MOD
        , hlCHROMA_SOURCE_MOD
        , hlCHROMA_FEYNHELL_TYPE_COUNT
        } hl_feynhell_type;
        for( u64 feynhell_type_index = 0; feynhell_type_index < hlCHROMA_FEYNHELL_TYPE_COUNT; ++feynhell_type_index ){
          hlXML_Result lam_real_query;
          hlXML_Result lam_imag_query;
          hlXML_Result lam_op_query  ;
          hlXML_Result lam_mom_query ;
          u64 opcode_offset;
          switch( feynhell_type_index ){
            case hlCHROMA_OLD_SINGLE_INVERSION_MOD:{
              //hlDEBUG_STR( "FH OLD\n" );
              lam_real_query = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_LAMBDA_REAL_PATH_OLD );
              lam_imag_query = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_LAMBDA_IMAG_PATH_OLD );
              lam_op_query   = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_OPERATOR_PATH_OLD    );  
              lam_mom_query  = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_MOMENTUM_PATH_OLD    );  
              opcode_offset = 0;
            }break;
            case hlCHROMA_INVERSION_MOD:{
              //hlDEBUG_STR( "FH NEW\n" );
              lam_real_query = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_LAMBDA_REAL_PATH );
              lam_imag_query = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_LAMBDA_IMAG_PATH );
              lam_op_query   = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_OPERATOR_PATH    );
              lam_mom_query  = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_MOMENTUM_PATH    );
              opcode_offset = 0;
            }break;
            case hlCHROMA_SOURCE_MOD:{
              //hlDEBUG_STR( "FH SOURCE\n" );
              lam_real_query = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_SOURCE_MOD_LAMBDA_REAL_PATH );
              lam_imag_query = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_SOURCE_MOD_LAMBDA_IMAG_PATH );
              lam_op_query   = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_SOURCE_MOD_OPERATOR_PATH    );
              lam_mom_query  = hlXMLQuery( &tmp_pool, &mem_pool, prop_xml, PROP_FEYNHELL_SOURCE_MOD_MOMENTUM_PATH    );
              opcode_offset = hlARRAYCOUNT( real_operator_identity_string ) - 1;
            }break;
          }
          hlASSERT( lam_real_query.count == lam_imag_query.count );
          hlASSERT( lam_real_query.count == lam_op_query.count   );
          hlASSERT( lam_real_query.count == lam_mom_query.count  );
          if( lam_real_query.count ){
            u64 previous_lambda_count     = correlator_params.props[prop_index].lambda_count;
            u64 new_lambda_addition_count = lam_real_query.count;
            correlator_params.props[prop_index].lambda_count = previous_lambda_count + new_lambda_addition_count;
            //printf( "[DEBUG] %d, %d -> %d\n", previous_lambda_count, new_lambda_addition_count, correlator_params.props[prop_index].lambda_count );
            for( u64 lam_index = previous_lambda_count; lam_index < correlator_params.props[prop_index].lambda_count; ++lam_index ){
              //hlDEBUG_STR( "ADDING MOD\n" );
              correlator_params.props[prop_index].feynman_hellmann_params[lam_index].lambda[0]     = lam_real_query.values[lam_index]->flt_val;
              correlator_params.props[prop_index].feynman_hellmann_params[lam_index].lambda[1]     = lam_imag_query.values[lam_index]->flt_val;
              correlator_params.props[prop_index].feynman_hellmann_params[lam_index].chroma_opcode = lam_op_query.values[lam_index]->num_val + opcode_offset;
              for( u64 mom_index = 0; mom_index < hlARRAYCOUNT( correlator_params.props[prop_index].feynman_hellmann_params[lam_index].lattice_momentum_transfer ); ++mom_index ){
                correlator_params.props[prop_index].feynman_hellmann_params[lam_index].lattice_momentum_transfer[mom_index] = lam_mom_query.values[lam_index]->num_arr_val.values[mom_index];
              }
            }
          }
        }
      }
      hlASSERT( mom2_max < hlARRAYCOUNT(mom2_max_to_momentum_number_count) );
      u64 np = mom2_max_to_momentum_number_count[mom2_max];
      for( u64 hadron_index = 0; hadron_index < hadron_numbers_count; ++hadron_index ){
        u64 hadron_number = hadron_numbers[ hadron_index ];
        correlator_params.hadron_number = hadron_number;
        for( u64 momentum_number = 0; momentum_number < np; ++momentum_number ){
          i64 zero_mom[] = { 0, 0, 0 };
          i64 *momentum = momentum_number_to_momenta[mom2_max][momentum_number];
          i64 mom2 = momentum[0]*momentum[0] + momentum[1]*momentum[1] + momentum[2]*momentum[2];
          i64 mom_sum = hlAbsi64( momentum[0] ) + hlAbsi64( momentum[1] ) + hlAbsi64( momentum[2] );
#if hlZERO_MOM
          if( mom2 == 0 ){
#elif hlONE_MOM
          if( mom2 == 1 ){
#endif
#if hlSIMPLE_MOM
          if( mom_sum <= 2 ){
#endif
      
          //if( momentum[ 0 ] == 0 && momentum[ 1 ] == 0 && momentum[ 2 ] >= 0 ){
            //-- GET OR GENERATE THE OBJECT TO COVER THE CORRELATOR TYPE --//
            correlator_params.lattice_momentum[0] = momentum[0];
            correlator_params.lattice_momentum[1] = momentum[1];
            correlator_params.lattice_momentum[2] = momentum[2];
            Correlator_Set *correlator_set = GetOrCopyCorrelatorSet( &correlator_hashmap, correlator_sets, correlator_params );

            //-- EXTRACT OUT THE CORRELATOR --//
            if( is_baryon_bin ){
              hlASSERT( hadron_number < hadron_bin_header->data_size / ( np * 2 * sizeof( r64 ) * nt ) );
              hlASSERT( hadron_number < hadron_trev_bin_header->data_size / ( np * 2 * sizeof( r64 ) * nt ) );
            }else{
              //printf( "%lu < %lu, %lu\n", hadron_number, hadron_bin_header->data_size, np*sizeof(r64)*nt );
              hlASSERT( hadron_number < hadron_bin_header->data_size / ( np * 1 * sizeof( r64 ) * nt ) );
            }
            u64 chunk_doubles = 2 * nt;
            u64 chunk_start   = hadron_number * np * nt * 2 * sizeof(r64) + momentum_number * nt * 2 * sizeof(r64);
            r64 *chunk        = hlPUSHARRAY( &mem_pool, chunk_doubles, r64 );
            hlMemCopy( hadron_bin_header->data + chunk_start, chunk, chunk_doubles * sizeof( r64 ) );
            BSwap64Array( chunk     , chunk_doubles );
            r64 *chunk_trev = 0;
            if( is_baryon_bin ){
              chunk_trev = hlPUSHARRAY( &mem_pool, chunk_doubles, r64 );
              hlMemCopy( hadron_trev_bin_header->data + chunk_start, chunk_trev, chunk_doubles * sizeof( r64 ) );
              BSwap64Array( chunk_trev, chunk_doubles );
            }
            //printf( "[DEBUG] %e, %e, %e, %e\n", chunk[0], chunk_trev[0], fwd_prop_cor_xml[0].values[0]->flt_arr_val.values[0], fwd_prop_cor_xml[1].values[0]->flt_arr_val.values[0] );
            if( !hlIsnan( chunk[0] ) && !(chunk_trev && hlIsnan( chunk_trev[0]) )
             && fwd_prop_cor_xml[ 0 ].values[ 0 ]->flt_arr_val.values[ 0 ] < 1e32 
             && fwd_prop_cor_xml[ 1 ].values[ 0 ]->flt_arr_val.values[ 0 ] < 1e32 ){ //--TODO BETTER CHECK FOR GOOD CORRELATOR --//
              //PrintSinkPair( sink_pair );
              //PrintSinkPair( collection->sink_pair );
              //Correlator *correlator = 0;
              //if( !correlator_set->first_correlator ){
              //  correlator = hlPUSHTYPE( &mem_pool, Correlator );
              //  correlator_set->first_correlator = correlator;
              //  correlator_set->last_correlator = correlator;
              //}else{
              //  correlator = hlPUSHTYPE( &mem_pool, Correlator );
              //  Correlator *end_correlator = correlator_set->last_correlator;
              //  end_correlator->next            = correlator;
              //  correlator_set->last_correlator = correlator;
              //}

              correlator.ident_string = ident_str;
              correlator.r64_count    = chunk_doubles;
              correlator.data         = chunk;
              correlator.data_trev     = chunk_trev;
              for( u64 prop_index = 0; prop_index < correlator_set->prop_count; ++prop_index ){
                correlator.forward_test_correlator[ prop_index ] = fwd_prop_cor_xml[ prop_index ].values[ 0 ]->flt_arr_val.values;
              }
              ++correlator_set->correlator_count;
              //PrintCorrelatorSet( *correlator_set );

              Correlator *stored_cor = hlPUSHTYPE( &mem_pool, Correlator );
              hlMemCopy( &correlator, stored_cor, sizeof(correlator) );
              if( !correlator_set->first_correlator ){
                correlator_set->first_correlator = stored_cor;
                correlator_set->last_correlator  = stored_cor;
              }else{
                Correlator *end_cor             = correlator_set->last_correlator;
                end_cor->next                   = stored_cor;
                correlator_set->last_correlator = stored_cor;
              }
            }else{
              if( first_bad_prop ){
                printf( "FOUND BAD CONFIGURATION PROP PAIR IN: %s\n", path );
                first_bad_prop = 0;
              }
              ++correlator_set->bad_correlators;
 
            }
          //}
#if hlZERO_MOM
          }
#elif hlONE_MOM
          }
#endif
#if hlSIMPLE_MOM
          }
#endif
        }
      }
    }
  }
  //for( u64 correlator_index = 0; correlator_index < NUM_CORRELATOR_SETS; ++correlator_index ){
  //  Correlator_Set set = correlator_sets[correlator_index];
  //  if( set.hash ){
  //    if( set.props[0].lambda_count && set.props[0].feynman_hellmann_params[0].lambda[0] > 0.0 && !set.props[1].lambda_count
  //     && set.props[0].sync_smear == hlSHELL_SMEARING && set.props[0].src_smear == hlSHELL_SMEARING && set.props[1].sync_smear == hlSHELL_SMEARING && set.props[1].src_smear == hlSHELL_SMEARING ){ 
  //      PrintCorrelatorSet( set );
  //    }
  //  }
  //}

  printf( "CONSTRUCTING JSON\n" );
  Correlator_Set *set = correlator_sets;
  for( Correlator_Set *set = correlator_sets; set->hash; ++set ){
    printf( "Found Next Collection ... " );
    printf( "%lu cors (bad:%lu)\n", set->correlator_count, set->bad_correlators );

    hlJSON root_object = {};
    hlJSONGenerateObject( &mem_pool, &root_object, 10 );

    if( set->hadron_type == hlBARYON ){
      //hlDEBUG_STR( "IT IS A BARYON\n" );
      hlJSONPushStringField( &root_object, "hadron", baryon_identity_string[ set->hadron_number ] );
      hlJSONPushStringField( &root_object, "spin"  , spin_identity_string  [ set->hadron_number ] );
    }else if( set->hadron_type == hlMESON ){
      //hlDEBUG_STR( "IT IS A MESON\n" );
      c8 *hadron_ident_1 = real_operator_identity_string[ 15 ]; //set->hadron_number / hlARRAYCOUNT( real_operator_identity_string ) ];
      c8 *hadron_ident_2 = real_operator_identity_string[ 15 ]; //set->hadron_number % hlARRAYCOUNT( real_operator_identity_string ) ];
      hlJSONPushStringField( &root_object, "hadron",  hlStrCatC( &mem_pool, hadron_ident_1, hadron_ident_2 ) );
      hlJSONPushStringField( &root_object, "spin"  , "unpol"   ); //--TODO DOESN'T ACTUALLY HAVE SPIN --//
    }

    hlJSON *dim_object = hlJSONPushArrayField( &mem_pool, &root_object, "dim", 4 );
    hlJSONPushNumberValue( dim_object, set->dim[0] );
    hlJSONPushNumberValue( dim_object, set->dim[1] );
    hlJSONPushNumberValue( dim_object, set->dim[2] );
    hlJSONPushNumberValue( dim_object, set->dim[3] );

    hlJSON *kappa_object = hlJSONPushArrayField( &mem_pool, &root_object, "kappas", 2 );
    hlJSONPushDoubleValue( kappa_object, set->sea_kappa1 );
    hlJSONPushDoubleValue( kappa_object, set->sea_kappa2 );

    hlJSON *momentum_object = hlJSONPushArrayField( &mem_pool, &root_object, "momentum", 3 );
    hlJSONPushNumberValue( momentum_object, set->lattice_momentum[ 0 ] ); 
    hlJSONPushNumberValue( momentum_object, set->lattice_momentum[ 1 ] ); 
    hlJSONPushNumberValue( momentum_object, set->lattice_momentum[ 2 ] ); 

    hlJSON *feyn_param_object = hlJSONPushArrayField( &mem_pool, &root_object, "propagators", set->prop_count );
    for( u64 prop_index = 0; prop_index < set->prop_count; ++prop_index ){
      newPropagator *prop = set->props + prop_index;
      hlJSON *prop_object = hlJSONPushObjectValue( &mem_pool, feyn_param_object, 4 );
      hlJSONPushDoubleField( prop_object, "kappa", prop->prop_kappa );
      
      u64 lambda_count = prop->lambda_count;
      if( lambda_count > 0 ){
        hlJSON *lambda_object = hlJSONPushArrayField( &mem_pool, prop_object, "feynman_hellmann", lambda_count );
        for( u64 lambda_index = 0; lambda_index < lambda_count; ++lambda_index ){
          Feynman_Hellmann_Params *param = prop->feynman_hellmann_params + lambda_index;
          hlJSON *feynman_hellmann_object = hlJSONPushObjectValue( &mem_pool, lambda_object, 3 );
          hlJSONPushDoubleField( feynman_hellmann_object, "lambda", param->lambda[1] != 0.0 ? param->lambda[1] : param->lambda[0] );
          hlJSONPushStringField( feynman_hellmann_object, "op"    , param->lambda[1] != 0.0 ? imag_operator_identity_string[param->chroma_opcode] : real_operator_identity_string[param->chroma_opcode] );
          hlJSON *momentum_transfer_object = hlJSONPushArrayField( &mem_pool, feynman_hellmann_object, "op_q", hlARRAYCOUNT( param->lattice_momentum_transfer ) );
          hlJSONPushNumberValue( momentum_transfer_object, param->lattice_momentum_transfer[ 0 ] );
          hlJSONPushNumberValue( momentum_transfer_object, param->lattice_momentum_transfer[ 1 ] );
          hlJSONPushNumberValue( momentum_transfer_object, param->lattice_momentum_transfer[ 2 ] );
        }
      }
    }

    hlJSON *extraction_object = hlJSONPushObjectField( &mem_pool, &root_object, "extraction", 2 );
    hlJSONPushNumberField( extraction_object, "bad_chunks" , set->bad_correlators  );
    hlJSONPushNumberField( extraction_object, "chunk_count", set->correlator_count );

    hlJSON *value_object = hlJSONPushArrayField( &mem_pool, &root_object, "keys", set->correlator_count );
    hlJSON *test_object  = hlJSONPushArrayField( &mem_pool, &root_object, "test", set->correlator_count );
    Correlator *current_cor = set->first_correlator;
    hlDEBUG_STR( "EXTRACTING CORS\n" );
    while( current_cor ){
      c8 *ident_str = hlPUSHARRAY( &mem_pool, 256, c8 );
      snprintf( ident_str, 256, "%s.%lu.%lu.%lu.%lu", current_cor->ident_string, current_cor->source[ 0 ]
              , current_cor->source[ 1 ], current_cor->source[ 2 ] , current_cor->source[ 3 ] );
      //printf( "%s.%lu.%lu.%lu.%lu\n", current_cor->ident_string, current_cor->source[ 0 ]
      //        , current_cor->source[ 1 ], current_cor->source[ 2 ] , current_cor->source[ 3 ] );
      hlJSONPushStringValue( value_object, ident_str ); hlJSONPushDoubleValue( test_object, current_cor->data[ 8 ] );
      //for( u64 double_index = 0; double_index < current_chunk->r64_count; ++double_index ){
      //  hlJSONPushDoubleValue( chunk_object, current_chunk->data[ double_index ] );
      //}
      current_cor = current_cor->next;
    }

    c8 lattice_directory[256];
    c8 run_identification_directory[256];
    c8 directories[2*256];
    c8 data_point_filename[256];
    c8 out_json_filename[3*256];
    c8 out_data_filename[3*256];

    c8 *hadron_ident_1 = 0;
    c8 *hadron_ident_2 = 0;
    b32 include_prop = 1;
    if( set->hadron_type == hlBARYON ){
      hadron_ident_1 = baryon_identity_string[ set->hadron_number ];
      hadron_ident_2 = spin_identity_string  [ set->hadron_number ];
    }else if( set->hadron_type == hlMESON ){
      hadron_ident_1 = real_operator_identity_string[ 15 ]; //set->hadron_number / hlARRAYCOUNT( real_operator_identity_string ) ];
      hadron_ident_2 = real_operator_identity_string[ 15 ]; //set->hadron_number % hlARRAYCOUNT( real_operator_identity_string ) ];
    }
    u64 folder1_index = snprintf( lattice_directory, sizeof( lattice_directory ), "%lu_%lu_%lu_%lu_%f_%f"
                                , set->dim[0], set->dim[1]
                                , set->dim[2], set->dim[3]
                                , set->sea_kappa1, set->sea_kappa2 );
    u64 folder2_index = 0;
    for( u64 prop_index = 0; prop_index < set->prop_count; ++prop_index ){
      newPropagator *prop = set->props + prop_index;
      if( set->sea_kappa1 == set->sea_kappa2 && hlAbs64( set->sea_kappa1 - prop->prop_kappa ) < R32_EPSILON ){
        folder2_index += snprintf( run_identification_directory + folder2_index, sizeof( run_identification_directory ) - folder2_index
                                 , "_prop%s%s", hlSmearingStrs[prop->sync_smear], hlSmearingStrs[prop->src_smear] );
      }else{
        folder2_index += snprintf( run_identification_directory + folder2_index, sizeof( run_identification_directory ) - folder2_index
                                 , "_prop%f%s%s", prop->prop_kappa, hlSmearingStrs[prop->sync_smear], hlSmearingStrs[prop->src_smear] );
      }
      if( prop->sync_smear == hlPOINT_SMEARING || prop->src_smear == hlPOINT_SMEARING ){
        hlDEBUG_STR( "ENCOUNTERED SMEARED SYNC OR SRC\n" );
        include_prop = 0;
      }
      for( u64 lambda_index = 0; lambda_index < prop->lambda_count; ++lambda_index ){
        Feynman_Hellmann_Params *param = prop->feynman_hellmann_params + lambda_index;
        r64 lambda;
        c8 *op_name = 0;
        if( param->lambda[1] != 0.0 ){
          lambda = param->lambda[1];
          op_name = imag_operator_identity_string[ param->chroma_opcode ];
        }else if( param->lambda[0] != 0.0 ){
          lambda = param->lambda[0];
          op_name = real_operator_identity_string[ param->chroma_opcode ];
        }
        if( op_name ){
          folder2_index += snprintf( run_identification_directory + folder2_index, sizeof( run_identification_directory ) - folder2_index
                                    , "_%s_%f_q%+li%+li%+li"
                                    , op_name, lambda
                                    , param->lattice_momentum_transfer[ 0 ], param->lattice_momentum_transfer[ 1 ], param->lattice_momentum_transfer[ 2 ] );
          hlASSERT( folder2_index < sizeof( run_identification_directory ) );
        }
     }
    }
    snprintf( data_point_filename, sizeof( data_point_filename ), "%s_%s_psync%+ld%+ld%+ld"
            , hadron_ident_1
            , hadron_ident_2
            , set->lattice_momentum[ 0 ]
            , set->lattice_momentum[ 1 ]
            , set->lattice_momentum[ 2 ] );
            
    snprintf( directories      , sizeof( directories       ), "%s/%s"        , lattice_directory, run_identification_directory                      );
    snprintf( out_data_filename, sizeof( out_data_filename ), "%s/%s/%s.dat" , lattice_directory, run_identification_directory, data_point_filename );
    snprintf( out_json_filename, sizeof( out_data_filename ), "%s/%s/%s.json", lattice_directory, run_identification_directory, data_point_filename );
    if( include_prop ){
      hlPrintProgramMemory( &tmp_pool, &mem_pool );
      mkdir( lattice_directory, 0777 );
      mkdir( directories, 0777 );
      printf( "Writting files: %s, %s\n", out_json_filename, out_data_filename );
      FILE *out_json_file = fopen( out_json_filename, "w+" );
      hlJSONPrint( out_json_file, &root_object );
      fclose( out_json_file );

      FILE *out_data_file = fopen( out_data_filename, "wb" );
      current_cor = set->first_correlator;
      while( current_cor ){
        //hlDEBUG_STR( "REV CHUNK\n" );
        fwrite( current_cor->data      , sizeof(r64), current_cor->r64_count, out_data_file );
        if( current_cor->data_trev ){
          //hlDEBUG_STR( "TREV CHUNK\n" );
          fwrite( current_cor->data_trev , sizeof(r64), current_cor->r64_count, out_data_file );
        }
        current_cor = current_cor->next;
      }    
      fclose( out_data_file );
      hlDEBUG_STR( "DATA FILE WRITTEN\n" );
    }
  }
    
  return( 0 );
}
