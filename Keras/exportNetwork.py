import argparse
from ROOT import *
from math import *
from keras.models import load_model
import datetime

def creat_ccmodel(model_file, outFile):
  #open file to store the model
  outputFile = open(outFile+'.h','w')

  model = load_model(model_file)
  
  outputFile.write("///Keras Model for c++ usage\n")
  outputFile.write("///Author: Miqueias Melo de Almeida\n\n")
  
  outputFile.write("///sigmoid output function\n")
  outputFile.write("double sigmoid_%s(double z){\n" % outFile)
  outputFile.write("\treturn 1./(1+TMath::Exp(-z));\n")
  outputFile.write("}\n\n")

  outputFile.write("///max function (ReLU)\n")
  outputFile.write("double relu_%s(double value){\n" % outFile)
  outputFile.write("\tdouble max_val = (value > 0)? value:0;\n")
  outputFile.write("\treturn max_val;\n")
  outputFile.write("}\n\n")

  outputFile.write("///PReLU function\n")
  outputFile.write("double prelu_%s(double slope, double value){\n" % outFile)
  outputFile.write("\tdouble max_val = (value > 0)? value:(value*slope);\n")
  outputFile.write("\treturn max_val;\n")
  outputFile.write("}\n\n")

  outputFile.write("///SeLU function\n")
  outputFile.write("double selu_%s(double value){\n" % outFile)
  outputFile.write("\tdouble alpha = 1.6732632423;\n")
  outputFile.write("\tdouble scale = 1.0507009873;\n")
  outputFile.write("\tdouble fvalue = (value >= 0)? (scale*value):(scale*alpha*(TMath::Exp(value)-1));\n")
  outputFile.write("\treturn fvalue;\n")
  outputFile.write("}\n\n")

  nlayers = len(model.layers)
  ilayer = nlayers
  layer_types = [-1 for i in range(nlayers)]
  #creates each neuron-input function
  outputFile.write("///Set of neuron functions with their respective weight and bias\n")
  while( ilayer-1 >=0 ):
    #g=model.layers[ilayer-1].get_config() #contains the weights and the bias
    h=model.layers[ilayer-1].get_weights() #contains the weights and the bias
    #print("Loading layer %i" % (ilayer-1))
    #print (g)
    #print (h)
    #h[0][i][j] --> weights from i-esimo input in the j-esimo neuron, h[1][j] --> bias for each neuron
  
    layer_type = len(h)
    layer_types[ilayer-1] = layer_type
    layer_ninputs = len(h[0])
    layer_nneurons = 1
    if(layer_type == 2):
      layer_nneurons = len(h[0][0])
    ineuron = layer_nneurons-1
    #print("Searching max absolute weight...")
  
    if(layer_type == 2):
      outputFile.write("double l%i_func_%s(std::vector<double> &inputs, int neuron){\n" % (ilayer,outFile))
      #outputFile.write("\tint n_inputs = inputs.size();\n")
    else:
      outputFile.write("double l%i_func_%s(int neuron){\n" % (ilayer,outFile))
    
    outputFile.write("\tdouble z = 0;\n")
    outputFile.write("\tswitch(neuron){\n")
    while( ineuron >=0 ): #look at each neuron (number of func/layer)
      if(layer_type == 2):
	outputFile.write("\t\tcase %i:\n" % (layer_nneurons-ineuron-1))
      n_inputs = layer_ninputs-1
      i_input = n_inputs
      while( i_input >=0 ):
	weight = 0
	if(layer_type == 2):
	  weight = h[0][n_inputs-i_input][layer_nneurons-ineuron-1]
	  outputFile.write("\t\t\tz += %.10f*inputs[%i];\n" %(weight,n_inputs-i_input))
	else:
	  weight = h[0][n_inputs-i_input]
	  outputFile.write("\t\tcase %i:\n" % (n_inputs-i_input))
	  outputFile.write("\t\t\tz = %.10f;\n" % weight)
	  outputFile.write("\t\tbreak;\n")
	i_input -= 1
      
      if(layer_type == 2):
	bias = h[1][layer_nneurons-ineuron-1]
	outputFile.write("\t\t\tz += %.10f;\n" % bias)      
	outputFile.write("\t\tbreak;\n")
      ineuron -= 1
    
    outputFile.write("\t}\n")
    outputFile.write("\n\treturn z;\n")
    outputFile.write("}\n")
    ilayer -= 1


  #create the network archictecture
  outputFile.write("\n\n")
  for ilayer in range(len(model.layers)):
    #g=layer.get_config()  #contains the network configuration
    h=model.layers[ilayer].get_weights() #contains the weights and the bias
    #print("Loading layer %i" % (ilayer))
    #print (g)
    #print (h)
    #h[0][i][j] --> weights from i-esimo input in the j-esimo neuron, h[1][j] --> bias for each neuron
  
    layer_type = len(h)
    layer_ninputs = len(h[0])
    layer_nneurons = 1
    if(layer_type == 2):
      layer_nneurons = len(h[0][0])
  
    #the first hidden layer
    if ilayer == 0:
      outputFile.write("double model_%s(std::vector<double> &inputs){\n" % outFile)
      outputFile.write("\tstd::vector<double> l%i_neurons_z_%s;\n" % (ilayer+1,outFile))
      outputFile.write("\tfor(int l%i_n=0; l%i_n<%i; ++l%i_n)\n" % (ilayer+1,ilayer+1,layer_nneurons,ilayer+1))
      if(layer_types[ilayer+1] == layer_type):
	#outputFile.write("\t\tl%i_neurons_z_%s.push_back( relu_%s(l%i_func_%s(inputs,l%i_n)) );\n" % (ilayer+1,outFile,outFile,ilayer+1,outFile,ilayer+1))
	outputFile.write("\t\tl%i_neurons_z_%s.push_back( selu_%s(l%i_func_%s(inputs,l%i_n)) );\n" % (ilayer+1,outFile,outFile,ilayer+1,outFile,ilayer+1))
      else:
	outputFile.write("\t\tl%i_neurons_z_%s.push_back( l%i_func_%s(inputs,l%i_n) );\n" % (ilayer+1,outFile,ilayer+1,outFile,ilayer+1))
      outputFile.write("//------ end layer %i ------\n\n" % (ilayer+1))
  
  
    #here's the complication (where out layer neurons appear)
    if ilayer > 0 and ilayer < nlayers-1:
      outputFile.write("\tstd::vector<double> l%i_neurons_z_%s;\n" % (ilayer+1,outFile))
      if(layer_types[ilayer+1] == layer_type):	
	outputFile.write("\tfor(int l%i_n=0; l%i_n<%i; ++l%i_n)\n" % (ilayer+1,ilayer+1,layer_nneurons,ilayer+1))
	#outputFile.write("\t\tl%i_neurons_z_%s.push_back( relu_%s(l%i_func_%s(l%i_neurons_z_%s,l%i_n)) );\n" % (ilayer+1,outFile,outFile,ilayer+1,outFile,ilayer,outFile,ilayer+1))
	outputFile.write("\t\tl%i_neurons_z_%s.push_back( selu_%s(l%i_func_%s(l%i_neurons_z_%s,l%i_n)) );\n" % (ilayer+1,outFile,outFile,ilayer+1,outFile,ilayer,outFile,ilayer+1))
      elif(layer_type == 1 and layer_types[ilayer+1] != layer_type):
	outputFile.write("\tfor(int l%i_n=0; l%i_n<%i; ++l%i_n)\n" % (ilayer+1,ilayer+1,layer_ninputs,ilayer+1))
	outputFile.write("\t\tl%i_neurons_z_%s.push_back( prelu_%s(l%i_func_%s(l%i_n),l%i_neurons_z_%s[l%i_n]) );\n" % (ilayer+1,outFile,outFile,ilayer+1,outFile,ilayer+1,ilayer,a,b,c,ilayer+1))
      else:
	outputFile.write("\tfor(int l%i_n=0; l%i_n<%i; ++l%i_n)\n" % (ilayer+1,ilayer+1,layer_ninputs,ilayer+1))
	outputFile.write("\t\tl%i_neurons_z_%s.push_back( l%i_func_%s(l%i_neurons_z_%s,l%i_n) );\n" % (ilayer+1,outFile,ilayer+1,outFile,ilayer,outFile,ilayer+1))
      outputFile.write("//------ end layer %i ------\n\n" % (ilayer+1))
  
  
    #last layer (output)
    if ilayer == nlayers-1:
      outputFile.write("\n\treturn sigmoid_%s( l%i_func_%s(l%i_neurons_z_%s,0) );\n" % (outFile,ilayer+1,outFile,ilayer,outFile))
      outputFile.write("}\n")
    
  outputFile.close()



def main(options):
  creat_ccmodel(options.infile, options.outfile)
  
  
  
if __name__ == '__main__':

 # Setup argument parser
 parser = argparse.ArgumentParser()

 # Add more arguments
 parser.add_argument("--infile", help="Name of txt input file with addresses of root files")
 parser.add_argument("--outfile", help="Name for model output file")

 # Parse default arguments
 options = parser.parse_args()
 main(options)
 
 print '!! Note: the model name is %s' % options.outfile
