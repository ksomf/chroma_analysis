//-------------------------------------------------------------------------------------------//
//-- Kim Somfleth <kim.somfleth@internode.on.net>                                          --//
//--                        LIME FILE COMMON USAGE FILE FOR READING THEM                   --//
//-------------------------------------------------------------------------------------------//

#if !defined( _CHROMA_LIME_READER_H )

#include "../hl.h"

typedef struct {
  u32 magic_number;
  u16 file_version;
  u16 magic_bits;
  u64 data_size;
  c8  type[ 128 ];
  u8  data[ 0 ];
} Lime_Header;

static void PrintLimeHeader( Lime_Header *lime_header ){
  printf( "M:%u, V:%hu, M:%hu, S:%lu, T:%.128s\n", lime_header->magic_number, lime_header->file_version, lime_header->magic_bits, lime_header->data_size, lime_header->type );
}

static Lime_Header *ReadLimeHeader( u8 **pos ){
  Lime_Header *result               = (Lime_Header *)*pos;
               result->magic_number = bswap32( result->magic_number );
               result->file_version = bswap16( result->file_version );
               result->magic_bits   = bswap16( result->magic_bits   );
               result->data_size    = bswap64( result->data_size    );
  //PrintLimeHeader( result );
  hlASSERT( result->magic_number == 1164413355 );
  *pos = ALIGN_U8_PTR_TO_8_BYTES( result->data + result->data_size );
  return( result );
}

static Lime_Header *ReadLimeHeader2( u8 **pos, c8 *debug ){
  printf( "DEBUG LIME READ: %s\n", debug );
  Lime_Header *result               = (Lime_Header *)*pos;
               result->magic_number = bswap32( result->magic_number );
               result->file_version = bswap16( result->file_version );
               result->magic_bits   = bswap16( result->magic_bits   );
               result->data_size    = bswap64( result->data_size    );
  //PrintLimeHeader( result );
  if( result->magic_number != 1164413355 ){
    printf( "BAD FILE: %s\n", debug );
  }
  hlASSERT( result->magic_number == 1164413355 );
  *pos = ALIGN_U8_PTR_TO_8_BYTES( result->data + result->data_size );
  return( result );
}

static c8 *GetConfigurationName( hlMem_Pool *tmp_pool, hlMem_Pool *mem_pool, hlXML *dir_xml, c8 *path ){
  c8* result = 0;
  //hlDEBUG_STR( "CHECKING FOR UNIT CONFIG\n" );
  if( hlXMLQuery( tmp_pool, tmp_pool, dir_xml, "Config_info/unit" ).count == 1 ){
    result = hlStrCopyC( mem_pool, "unit" );
  }else{
    //hlDEBUG_STR( "UNIT CONFIG NOT FOUND\n" );
    hlXML_Result ildg_config_name_xml        = hlXMLQuery( tmp_pool, tmp_pool, dir_xml, "Config_info/ILDG/LFN" ); 
    hlXML_Result szinqio_inv_config_name_xml = hlXMLQuery( tmp_pool, tmp_pool, dir_xml, "Input/lime_file" );
    hlXML_Result szinqio_gen_config_name_xml = hlXMLQuery( tmp_pool, tmp_pool, dir_xml, "Input/chroma/Cfg/cfg_file" );
    c8 *query_result;
    //printf( "%lu, %lu, %lu\n", ildg_config_name_xml.count, szinqio_inv_config_name_xml.count, szinqio_gen_config_name_xml.count );
    if( ildg_config_name_xml.count == 1 || szinqio_inv_config_name_xml.count == 1 ){
      //printf( "ILDG CONFIG OR SZINQIO HADRON SPECTROSCOPY OUTPUT\n" );
      hlASSERT( (szinqio_inv_config_name_xml.count == 1 && szinqio_inv_config_name_xml .values[ 0 ]->type == hlXML_STR )
            || (ildg_config_name_xml.count == 1 && ildg_config_name_xml.values[ 0 ]->type == hlXML_STR) );
      query_result = ildg_config_name_xml.count ? ildg_config_name_xml.values[ 0 ]->str 
                                                : szinqio_inv_config_name_xml.values[ 0 ]->str;
      query_result = hlStrAfterLastC( '/', query_result );
      //printf( "GENERATING CONFIG NAME FROM: %s\n", query_result );
      c8 *before_gen    = hlStrAfterC( '.', query_result );
      c8 *config_start  = hlStrAfterC( '.', before_gen );
      c8 *gen_end       = config_start;
      if( hlStrStartsWithC( "b5", config_start ) ){ //REMOVE OCCURENSES OF BETA AND LIKE FROM FILE
        config_start = hlStrAfterC( '.', config_start );
      }
      c8 tokken_for_next = '_';
      if( !hlStrCountTokkenC( config_start, tokken_for_next ) ){
        tokken_for_next = '.';
      }
      c8 *after_config  = hlStrAfterC( tokken_for_next, config_start );
      c8 *end_config    = after_config - 1;
      
      //printf( "%lu\n", (u64)(end_config-config_start) );
      c8 *config = hlStrCopy( tmp_pool, (u64)(end_config-config_start), config_start );
      u64 config_u64 = strtoul( config, 0, 10 );
      //printf( "%lu ( %s -> %lu )\n", (u64)(end_config-config_start), config, config_u64 );
      result = hlStrCopyC( mem_pool, query_result );
      snprintf( result + (u64)(gen_end - query_result), 6, "%05lu", config_u64 );
      
      //printf( "%s\n", result );
      //hlASSERT( !"STOP FOR NOW" );
    }else if( szinqio_gen_config_name_xml.count == 1 ){
      hlASSERT( szinqio_gen_config_name_xml.values[ 0 ]->type == hlXML_STR );
      //printf( "SZINQIO CONFIG GEN OUTPUT\n" );
      query_result = szinqio_gen_config_name_xml.values[ 0 ]->str;
      query_result = hlStrAfterLastC( '/', query_result );
      c8 *before_suffix = hlStrAfterLastC( '.', query_result ) - 1;
      c8 *number_start  = hlStrAfterLastC( 'e', before_suffix );
      u64 number = strtoul( number_start, 0, 10 );
      
      result = hlStrCopy( mem_pool, (u64)(before_suffix-query_result) + 9, query_result );
      snprintf( result + (u64)(before_suffix-query_result), 3+6, ".1.%05lu", number );
      //printf( "%s <- %s\n", result, query_result );
      
    }else if( path ){
      query_result = path;
      query_result = hlStrAfterLastC( '/', query_result );
      //printf( "GENERATING CONFIG NAME FROM: %s\n", query_result );
      c8 *before_gen    = hlStrAfterC( '.', query_result );
      c8 *config_start  = hlStrAfterC( '.', before_gen );
      c8 *after_config  = hlStrAfterC( '_', config_start );
      c8 *end_config    = after_config - 1;
      
      //printf( "%lu\n", (u64)(end_config-config_start) );
      c8 *config = hlStrCopy( tmp_pool, (u64)(end_config-config_start), config_start );
      u64 config_u64 = strtoul( config, 0, 10 );
      //printf( "%lu ( %s -> %lu )\n", (u64)(end_config-config_start), config, config_u64 );
      result = hlStrCopyC( mem_pool, query_result );
      snprintf( result + (u64)(config_start - query_result), 6, "%05lu", config_u64 );
      
      //printf( "%s\n", result );
    }else{
      hlASSERT( !"NOT ABLE TO FIND CONFIGURATION NAME" );
    }
  }
  
  return( result );
}

#define _CHROMA_LIME_READER_H
#endif
