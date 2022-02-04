<?xml version="1.0"?>
<chroma>
<annotation>
  Hadron Spectrum Input
</annotation>
<Param>
  <InlineMeasurements>
{{#wilson_flows}}
<elem>
  <annotation>WilsonFlow the gauge field</annotation>
  <Name>WILSON_FLOW</Name>
  <Frequency>1</Frequency>
  <Param>
    <version>1</version>
    <nstep>{{incwsteps}}</nstep>
    <wtime>{{incwtime}}</wtime>
    <t_dir>3</t_dir>
  </Param>
  <NamedObject>
    <gauge_in>{{#start}}{{> gauge_id}}{{/start}}{{^start}}default_gauge_field{{/start}}</gauge_in>
    <gauge_out>{{> gauge_id}}</gauge_out>
  </NamedObject>
</elem>

<elem>
  <Name>QIO_WRITE_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> gauge_id}}</object_id>
    <object_type>Multi1dLatticeColorMatrix</object_type>
  </NamedObject>
  <File>
    <file_name>{{{TMP_DIR}}}/{{{TRAJ}}}_{{> gauge_id}}.lime</file_name>
    <file_volfmt>SINGLEFILE</file_volfmt>
  </File>
</elem>
{{#start}}<elem>
  <Name>ERASE_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> gauge_id}}</object_id>
  </NamedObject>
</elem>{{/start}}
{{/wilson_flows}}

  </InlineMeasurements>
   <nrow>24 24 24 48</nrow>
</Param>

<RNG>
  <Seed>
    <elem>11</elem>
    <elem>11</elem>
    <elem>11</elem>
    <elem>0</elem>
  </Seed>
</RNG>

<Cfg>
    <cfg_type>{{cfg_type}}</cfg_type>
    <cfg_file>{{{CONF}}}</cfg_file>
</Cfg>
</chroma>
