<?xml version="1.0"?>
<chroma>
<annotation>
  Hadron Spectrum Input
</annotation>
<Param> 
  <InlineMeasurements>

{{#readins}}
<elem>
  <Name>QIO_READ_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> prop_id}}_smeared_sink</object_id>
    <object_type>LatticePropagator</object_type>
  </NamedObject>
  <File>
    <file_name>{{{TMP_DIR}}}/{{> prop_id}}_smeared_sink.lime</file_name>
  </File>
</elem>
{{/readins}}

<elem>
  <Name>BARYON_SPECTRUM-QCDSF</Name>
  <Frequency>1</Frequency>
  <Param>
    <version>1</version>
    <fwdbwd_average>false</fwdbwd_average>
    <time_rev>true</time_rev>
    <mom2_max>{{mom2_max}}{{^mom2_max}}0{{/mom2_max}}</mom2_max>
    <avg_equiv_mom>false</avg_equiv_mom>
    <xml>false</xml>
    <lime>true</lime>
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <sink_pairs>

{{#combinations}}
<elem>
  <first_id>{{#prop1}}{{> prop_id}}{{/prop1}}_smeared_sink</first_id>
  <second_id>{{#prop2}}{{> prop_id}}{{/prop2}}_smeared_sink</second_id>
</elem>
{{/combinations}}
</sink_pairs>
</NamedObject>
<lime_file>{{{TMP_DIR}}}/baryon_{{{TRAJ}}}_{{#slice}}{{slice}}{{/slice}}{{^slice}}0.000{{/slice}}.lime</lime_file>
<xml_file>{{{TMP_DIR}}}/baryon_{{{TRAJ}}}_{{#slice}}{{slice}}{{/slice}}{{^slice}}0.000{{/slice}}.xml</xml_file>
</elem>

<elem>
  <Name>MESON_SPECTRUM-QCDSF{{#SMALL_MESON}}-SMALL{{/SMALL_MESON}}</Name>
  <Frequency>1</Frequency>
  <Param>
    <version>1</version>
    <fwdbwd_average>false</fwdbwd_average>
    <time_rev>true</time_rev>
    <mom2_max>{{mom2_max}}{{^mom2_max}}0{{/mom2_max}}</mom2_max>
    <avg_equiv_mom>false</avg_equiv_mom>
    <xml>false</xml>
    <lime>true</lime>
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <sink_pairs>

{{#combinations}}
<elem>
  <first_id>{{#prop1}}{{> prop_id}}{{/prop1}}_smeared_sink</first_id>
  <second_id>{{#prop2}}{{> prop_id}}{{/prop2}}_smeared_sink</second_id>
</elem>
{{/combinations}}

</sink_pairs>
</NamedObject>
<lime_file>{{{TMP_DIR}}}/meson_{{{TRAJ}}}_{{#slice}}{{slice}}{{/slice}}{{^slice}}0.000{{/slice}}.lime</lime_file>
<xml_file>{{{TMP_DIR}}}/meson_{{{TRAJ}}}_{{#slice}}{{slice}}{{/slice}}{{^slice}}0.000{{/slice}}.xml</xml_file>
</elem>

  </InlineMeasurements>
   <nrow>{{{NX}}} {{{NY}}} {{{NZ}}} {{{NT}}}</nrow>
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

