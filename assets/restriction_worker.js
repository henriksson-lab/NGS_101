'use strict';
importScripts('sequence_annotation.js','restriction_analysis.js');
self.onmessage = function(event) {
  try {
    const {sequence,enzymes,circular}=event.data;
    const result=self.RestrictionAnalysis.analyze(sequence,enzymes,{circular,progress:(done,total)=>self.postMessage({type:'progress',done,total})});
    self.postMessage({type:'result',result});
  } catch(error) { self.postMessage({type:'error',error:error.message || String(error)}); }
};
