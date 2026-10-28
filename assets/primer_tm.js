/* Published DNA/DNA thermodynamics, independent of vendor calculator code.
 * SantaLucia 1998 table 2 / eq. 3; Owczarzy 2008 eqs. 16, 18–20 and decision tree.
 * Input concentrations: nM oligos, mM ions and TOTAL dNTPs. Internals use mol/L.
 */
(function(root,factory){
  const api=factory(typeof module==='object'&&module.exports?require('./sequence_annotation.js'):root.SequenceAnnotation);
  if(typeof module==='object'&&module.exports)module.exports=api;else root.PrimerTm=api;
}(typeof globalThis!=='undefined'?globalThis:this,function(sequenceTools){
  'use strict';
  function finite(value,name,min,max){if(typeof value!=='number'||!Number.isFinite(value)||value<min||value>max)throw new Error(`Enter ${name} between ${min} and ${max}.`);return value;}
  function normalize(text){const s=sequenceTools.normalizeInput(text);if(!/^[ACGT]{2,200}$/.test(s))throw new Error('Use 2–200 unmodified A/C/G/T bases per primer.');return s;}
  function freeMagnesium(mg_mM,dntpTotal_mM){
    finite(mg_mM,'Mg²⁺ (mM)',0,100);finite(dntpTotal_mM,'total dNTPs (mM)',0,100);
    // One-to-one Mg–dNTP equilibrium; association constant 3 × 10^4 M^-1.
    const mg=mg_mM/1000,d=dntpTotal_mM/1000,k=3e4,b=1+k*(d-mg),root=Math.sqrt(b*b+4*k*mg);
    return (b>=0?2*mg/(b+root):(root-b)/(2*k))*1000;
  }
  function saltCorrection(sequence,monovalent_mM,mg_mM,dntpTotal_mM){
    finite(monovalent_mM,'monovalent ions (mM)',0,1000);
    const mon=monovalent_mM/1000,mg=freeMagnesium(mg_mM,dntpTotal_mM)/1000,gc=[...sequence].filter(b=>b==='G'||b==='C').length/sequence.length;
    if(mon===0&&mg===0)throw new Error('The model needs monovalent ions or magnesium.');
    const ratio=mon===0?Infinity:Math.sqrt(mg)/mon;
    let correction,regime;
    if(ratio<.22){
      const log=Math.log(mon);correction=(4.29*gc-3.95)*1e-5*log+9.4e-6*log*log;regime='monovalent';
    }else{
      let a=3.92e-5,d=1.42e-5,g=8.31e-5;
      if(ratio<6){const log=Math.log(mon);a*=.843-.352*Math.sqrt(mon)*log;d*=1.279-.00403*log-.00803*log*log;g*=.486-.258*log+.00525*log**3;regime='mixed';}
      else regime='magnesium';
      const log=Math.log(mg);correction=a-9.11e-6*log+gc*(6.26e-5+d*log)+(-4.82e-4+5.25e-4*log+g*log*log)/(2*(sequence.length-1));
    }
    return {inverseKelvin:correction,regime,freeMg_mM:mg*1000,ratio};
  }
  function correctTemperature(referenceKelvin,sequence,conditions){
    finite(referenceKelvin,'reference temperature (K)',1,1000);
    const salt=saltCorrection(sequence,conditions.monovalent_mM,conditions.mg_mM,conditions.dntpTotal_mM);
    const kelvin=1/(1/referenceKelvin+salt.inverseKelvin);
    if(!Number.isFinite(kelvin)||kelvin<=0)throw new Error('These conditions are outside the usable temperature model.');
    return {...salt,tm_C:kelvin-273.15};
  }
  function calculate(text,conditions,model,tail=0){
    const fullSequence=normalize(text);if(!Number.isInteger(tail)||tail<0||tail>fullSequence.length-2)throw new Error('The 5′ tail must leave at least two binding bases.');
    const sequence=fullSequence.slice(tail),primer=finite(conditions.primer_nM,'primer concentration (nM)',.001,1e6),target=finite(conditions.target_nM,'target concentration (nM)',0,primer);
    if(!model||!model.nearest_neighbours)throw new Error('Missing thermodynamic parameters.');
    let dh=0,ds=0;
    for(let i=0;i<sequence.length-1;i++){const pair=model.nearest_neighbours[sequence.slice(i,i+2)];dh+=pair[0];ds+=pair[1];}
    for(const base of [sequence[0],sequence[sequence.length-1]]){const [h,s]=base==='G'||base==='C'?model.terminal_gc:model.terminal_at;dh+=h;ds+=s;}
    const selfComplementary=sequence===sequenceTools.revcomp(sequence);
    // SantaLucia eq. 3: C_A - C_B/2; self-complementary strands use total C_T.
    const effective_nM=selfComplementary?primer+target:primer-target/2;
    if(selfComplementary)ds+=model.symmetry_entropy;
    const referenceKelvin=dh*1000/(ds+model.gas_constant*Math.log(effective_nM*1e-9));
    const corrected=correctTemperature(referenceKelvin,sequence,conditions),warnings=[];
    if(sequence.length<10||sequence.length>60)warnings.push('Binding length is outside the 10–60 bp range discussed for the salt model.');
    if(corrected.regime!=='monovalent'&&(corrected.freeMg_mM<.5||corrected.freeMg_mM>600))warnings.push('Free magnesium is outside the published 0.5–600 mM magnesium-fit range.');
    return {...corrected,fullSequence,sequence,tail,length:sequence.length,gc_percent:100*[...sequence].filter(b=>b==='G'||b==='C').length/sequence.length,dh_kcal:dh,ds_cal:ds,referenceTm_C:referenceKelvin-273.15,effective_nM,selfComplementary,warnings};
  }
  return {normalize,freeMagnesium,saltCorrection,correctTemperature,calculate};
}));
