<?xml version="1.0"?>
<chroma>
<annotation>
  Hadron Spectrum Input
</annotation>
<Param> 
  <InlineMeasurements>
<elem>
  <Name>MAKE_SOURCE</Name>
  <Frequency>1</Frequency>
  <Param>
    <version>6</version>
    <Source>
      <version>3</version>
      <SourceType>SHELL_SOURCE</SourceType>
      <j_decay>3</j_decay>
      <t_srce>{{{SX}}} {{{SY}}} {{{SZ}}} {{{ST}}}</t_srce>
      <quark_smear_lastP>false</quark_smear_lastP>
      <SmearingParam>
        <wvf_kind>GAUGE_INV_JACOBI</wvf_kind>
        <wvf_param>0.21</wvf_param>
        <wvfIntPar>75</wvfIntPar>
        <no_smear_dir>3</no_smear_dir>
      </SmearingParam>
      <Displacement>
        <version>1</version>
        <DisplacementType>NONE</DisplacementType>
      </Displacement>
    </Source>
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <source_id>sh_source_1</source_id>
  </NamedObject>
</elem>

<!-- THE UNPERTUORPED PROPOGATOR -->
<elem>
  <Name>PROPAGATOR</Name>
  <Frequency>1</Frequency>
  <Param>
    <version>10</version>
    <quarkSpinType>FULL</quarkSpinType>
    <obsvP>false</obsvP>
    <numRetries>1</numRetries>
    <FermionAction>
      <FermAct>{{^dds}}{{fermion.action}}{{/dds}}{{#dds}}UNPRECONDITIONED_SLRC{{/dds}}</FermAct>
      <Kappa>{{kappa}}</Kappa>                    <!-- light mass -->
      <clovCoeff>{{fermion.csw}}</clovCoeff>{{^fermion.noslic}}
      <FermState>
        <Name>SLIC_FERM_STATE</Name>
        <n_smear>1</n_smear>
        <rho>0.1</rho>
        <orthog_dir>5</orthog_dir>{{/fermion.noslic}}
        <FermionBC>
          <FermBC>SIMPLE_FERMBC</FermBC>
          <boundary>1 1 1 -1</boundary>
        </FermionBC>{{^fermion.noslic}}
      </FermState>{{/fermion.noslic}}
    </FermionAction>
       {{#gpus}}<InvertParam>
         <invType>BICGSTAB_INVERTER</invType>
         <RsdBiCGStab>5.0e-11</RsdBiCGStab>
         <MaxBiCGStab>10000</MaxBiCGStab>
       </InvertParam>{{/gpus}}
         {{^gpus}}{{#dds}}<InvertParam>
            <invType>DDS_INVERTER</invType>
            <Kappa>{{kappa}}</Kappa>
            <Csw>{{fermion.csw}}</Csw>
            <Nmr>5</Nmr>  
            <Ncy>8</Ncy>
            <Nkv>18</Nkv>
            <DeflatedNV>4</DeflatedNV>
            <BlockSize>{{geometry}}</BlockSize>
            <RsdDDS>1.0e-12</RsdDDS>
            <MaxDDS>1000</MaxDDS>
            <BlkRel>0</BlkRel>
            <SLRC>1</SLRC>
          </InvertParam>{{/dds}}
         {{^dds}}{{#simple_invert}} <InvertParam>
           <invType>BICGSTAB_INVERTER</invType>
           <RsdBiCGStab>5.0e-11</RsdBiCGStab>
           <MaxBiCGStab>10000</MaxBiCGStab>
        </InvertParam>{{/simple_invert}}
         {{^simple_invert}}<InvertParam>
           <invType>RELIABLE_BICGSTAB_MP_SLRC_INVERTER-QCDSF</invType>
           <MaxIter>10000</MaxIter>
           <RsdTarget>1e-12</RsdTarget>
           <Delta>0.01</Delta>
           <CloverParams>
             <Kappa>{{kappa}}</Kappa>
             <clovCoeff>{{fermion.csw}}</clovCoeff>
             <FermState>
               <Name>SLIC_FERM_STATE</Name>
               <version>2</version>
               <n_smear>1</n_smear>
               <rho>0.1</rho>
               <smear_in_this_dirP>1 1 1 1</smear_in_this_dirP>
               <FermionBC>
                 <FermBC>SIMPLE_FERMBC</FermBC>
                 <boundary>1 1 1 -1</boundary>
               </FermionBC>
             </FermState>
           </CloverParams>
         </InvertParam>{{/simple_invert}}{{/dds}}{{/gpus}}
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <source_id>sh_source_1</source_id>
    <prop_id>{{> prop_id}}</prop_id>
  </NamedObject>
</elem>
<elem>
  <Name>SINK_SMEAR</Name>
  <Frequency>1</Frequency>
  <Param>
    <version>5</version>
    <Sink>
      <version>2</version>
      <SinkType>SHELL_SINK</SinkType>
      <j_decay>3</j_decay>
      <SmearingParam>
        <wvf_kind>GAUGE_INV_JACOBI</wvf_kind>
        <wvf_param>0.21</wvf_param>
        <wvfIntPar>75</wvfIntPar>
        <no_smear_dir>3</no_smear_dir>
      </SmearingParam>
      <Displacement>
        <version>1</version>
        <DisplacementType>NONE</DisplacementType>
      </Displacement>
    </Sink>
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <prop_id>{{> prop_id}}</prop_id>
    <smeared_prop_id>{{> prop_id}}_smeared_sink</smeared_prop_id>
  </NamedObject>
</elem>
<elem>
  <Name>QIO_WRITE_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> prop_id}}_smeared_sink</object_id>
    <object_type>LatticePropagator</object_type>
  </NamedObject>
  <File>
    <file_name>{{{TMP_DIR}}}/{{> prop_id}}_smeared_sink.lime</file_name>
    <file_volfmt>SINGLEFILE</file_volfmt>
  </File>
</elem>
<!-- END PREAMBLE -->

{{#propagators}}
<elem>
  <Name>PROPAGATOR</Name>
  <Frequency>1</Frequency>
  <Param>
    <version>10</version>
    <quarkSpinType>FULL</quarkSpinType>
    <obsvP>false</obsvP>
    <numRetries>1</numRetries>
    <FermionAction>
      <FermAct>UNPREC_SLRC_FEYNHELL</FermAct>
      <Kappa>{{kappa}}</Kappa>                    <!-- light mass -->
      <FeynHellParam>
        {{#fh-params}}
        <elem>
          <LambdaReal>{{lam_real}}</LambdaReal>
          <LambdaImag>{{lam_imag}}</LambdaImag>
          <Operator>{{op_num}}</Operator>
          <Momentum>{{#q}}{{.}} {{/q}}</Momentum>
          <Source>{{{SX}}} {{{SY}}} {{{SZ}}} {{{ST}}}</Source>
          <NoiseReal>1.0</NoiseReal>
          <NoiseImag>0.0</NoiseImag>
        </elem>
        {{/fh-params}}
      </FeynHellParam>
      <clovCoeff>{{fermion.csw}}</clovCoeff>
      <FermState>
        <Name>SLIC_FERM_STATE</Name>
        <n_smear>1</n_smear>
        <rho>0.1</rho>
        <orthog_dir>5</orthog_dir>
        <FermionBC>
          <FermBC>SIMPLE_FERMBC</FermBC>
          <boundary>1 1 1 -1</boundary>
        </FermionBC>
      </FermState>
    </FermionAction>
       {{#gpus}}<InvertParam>
           <invType>BICGSTAB_INVERTER</invType>
           <RsdBiCGStab>5.0e-11</RsdBiCGStab>
           <MaxBiCGStab>10000</MaxBiCGStab>
       </InvertParam>{{/gpus}}
       {{^gpus}} <InvertParam>
            <invType>RELIABLE_BICGSTAB_MP_SLRC_FEYNHELL_INVERTER-QCDSF</invType>
            <MaxIter>10000</MaxIter>
            <RsdTarget>1e-12</RsdTarget>
            <Delta>0.01</Delta>
            <CloverParams>
              <Kappa>{{kappa}}</Kappa>
              <clovCoeff>{{fermion.csw}}</clovCoeff>
              <FermState>
                <Name>SLIC_FERM_STATE</Name>
                <version>2</version>
                <n_smear>1</n_smear>
                <rho>0.1</rho>
                <smear_in_this_dirP>1 1 1 1</smear_in_this_dirP>
                <FermionBC>
                  <FermBC>SIMPLE_FERMBC</FermBC>
                  <boundary>1 1 1 -1</boundary>
                </FermionBC>
              </FermState>
            </CloverParams>
            <FeynHellParam>
              {{#fh-params}}
              <elem>
                <LambdaReal>{{lam_real}}</LambdaReal>
                <LambdaImag>{{lam_imag}}</LambdaImag>
                <Operator>{{op_num}}</Operator>
                <Momentum>{{#q}}{{.}} {{/q}}</Momentum>
                <Source>{{{SX}}} {{{SY}}} {{{SZ}}} {{{ST}}}</Source>
                <NoiseReal>1.0</NoiseReal>
                <NoiseImag>0.0</NoiseImag>
              </elem>
              {{/fh-params}}
            </FeynHellParam>
          </InvertParam>{{/gpus}}
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <source_id>sh_source_1</source_id>
    <prop_id>{{> prop_id}}</prop_id>
    {{#initial_guess}}<initial_guess_id>{{> prop_id}}</initial_guess_id>{{/initial_guess}}{{^initial_guess}}zero{{/initial_guess}}
{{#prior_perturbed_guess}}    <prior_perturbed_guess_id>{{> prop_id}}</prior_perturbed_guess_id>
    <prior_perturbed_factor>{{factor}}</prior_perturbed_factor>
{{/prior_perturbed_guess}}
  </NamedObject>
</elem>
<elem>
  <Name>SINK_SMEAR</Name>
  <Frequency>1</Frequency>
  <Param>
    <version>5</version>
    <Sink>
      <version>2</version>
      <SinkType>SHELL_SINK</SinkType>
      <j_decay>3</j_decay>
      <SmearingParam>
        <wvf_kind>GAUGE_INV_JACOBI</wvf_kind>
        <wvf_param>0.21</wvf_param>
        <wvfIntPar>75</wvfIntPar>
        <no_smear_dir>3</no_smear_dir>
      </SmearingParam>
      <Displacement>
        <version>1</version>
        <DisplacementType>NONE</DisplacementType>
      </Displacement>
    </Sink>
  </Param>
  <NamedObject>
    <gauge_id>default_gauge_field</gauge_id>
    <prop_id>{{> prop_id}}</prop_id>
    <smeared_prop_id>{{> prop_id}}_smeared_sink</smeared_prop_id>
  </NamedObject>
</elem>
<elem>
  <Name>QIO_WRITE_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> prop_id}}_smeared_sink</object_id>
    <object_type>LatticePropagator</object_type>
  </NamedObject>
  <File>
    <file_name>{{{TMP_DIR}}}/{{> prop_id}}_smeared_sink.lime</file_name>
    <file_volfmt>SINGLEFILE</file_volfmt>
  </File>
</elem>
{{#erase_objects}}
<elem>
  <Name>ERASE_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> prop_id}}</object_id>
  </NamedObject>
</elem>
{{/erase_objects}}
{{^unified_run}}
<elem>
  <Name>ERASE_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> prop_id}}_smeared_sink</object_id>
  </NamedObject>
</elem>
{{/unified_run}}
{{/propagators}}
{{#unified_run}}
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
{{/unified_run}}

<!-- ERASE THE UNPERTURBED PROP -->
<elem>
  <Name>ERASE_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> prop_id}}</object_id>
  </NamedObject>
</elem>
<elem>
  <Name>ERASE_NAMED_OBJECT</Name>
  <Frequency>1</Frequency>
  <NamedObject>
    <object_id>{{> prop_id}}_smeared_sink</object_id>
  </NamedObject>
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

