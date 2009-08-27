lib_LTLIBRARIES = librnnlib.la
bin_PROGRAMS = rnnlib
rnnlib_SOURCES = Main.cpp
rnnlib_LDADD = librnnlib.la
librnnlib_la_includedir = $(includedir)/rnnlib
librnnlib_la_include_HEADERS = ActivationFunctions.hpp  DecodingLayer.hpp      Mdrnn.hpp            SamplingLayer.hpp \
AutoregressionLayer.hpp  FullConnection.hpp     MultiArray.hpp       SeqBuffer.hpp \
BiasLayer.hpp            GatherLayer.hpp        MultilayerNet.hpp    SoftmaxLayer.hpp \
BlockLayer.hpp           GradientCheck.hpp      Named.hpp            SteepestDescent.hpp \
ClassificationLayer.hpp  Helpers.hpp            NetcdfDataset.hpp    StringAlignment.hpp \
CollapseLayer.hpp        InputLayer.hpp         NeuronLayer.hpp      Trainer.hpp \
ConfigFile.hpp           Layer.hpp              ObjectStore.hpp      TranscriptionLayer.hpp \
Connection.hpp           LinearOutputLayer.hpp  Optimiser.hpp        WeightContainer.hpp \
CopyConnection.hpp       list_of.hpp            OutputLayer.hpp      Word.hpp \
DataExporter.hpp         LogDouble.hpp          Random.hpp           WordTranscriptionLayer.hpp \
DataSequence.hpp         LstmLayer.hpp          RegressionLayer.hpp \
DatasetErrors.hpp        Matrix.hpp             Rprop.hpp
librnnlib_la_SOURCES = Random.cpp DataExporter.cpp
librnnlib_la_LDFLAGS = -version-info 2:0:0
