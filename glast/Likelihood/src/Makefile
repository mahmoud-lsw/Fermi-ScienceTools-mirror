HD_COMPONENT_NAME       = glast

HD_COMPONENT_VERS       =

HD_LIBRARY_ROOT         = Likelihood

HD_LIBRARY_SRC_cxx      = Accumulator.cxx \
                          AppHelpers.cxx \
                          BandFunction.cxx \
                          BinnedExposure.cxx \
                          BinnedHealpixExposure.cxx \
                          BinnedLikelihood.cxx \
                          BinnedLikelihood2.cxx \
                          BrokenPowerLaw2.cxx \
                          BrokenPowerLaw3.cxx \
                          BrokenPowerLawExpCutoff.cxx \
                          Composite2.cxx \
                          CompositeLikelihood.cxx \
                          Convolve.cxx \
                          CountsMap.cxx \
                          CountsSpectra.cxx \
                          DiffRespIntegrand.cxx \
                          DiffRespIntegrand2.cxx \
                          DiffRespNames.cxx \
                          DiffuseSource.cxx \
                          Drm.cxx \
                          EblAtten.cxx \
                          EnergyBand.cxx \
                          EquinoxRotation.cxx \
                          Event.cxx \
                          EventContainer.cxx \
                          ExpCutoff.cxx \
                          ExpCutoffSEDPeak.cxx \
                          ExposureCube.cxx \
                          ExposureMap.cxx \
                          FileFunction.cxx \
                          FluxBuilder.cxx \
                          HistND.cxx \
                          LikeExposure.cxx \
                          LogGaussian.cxx \
                          LogLike.cxx \
                          LogNormal.cxx \
                          LogNormalMuDist.cxx \
                          LogParabola.cxx \
                          MapBase.cxx \
                          MapCubeFunction2.cxx \
                          MeanPsf.cxx \
                          ModelMap.cxx \
                          MultipleBrokenPowerLaw.cxx \
                          Npred.cxx \
                          OneSourceFunc.cxx \
                          OptEM.cxx \
                          PiecewisePowerLaw.cxx \
                          Pixel.cxx \
                          PointSource.cxx \
                          PowerLaw2.cxx \
                          PowerLawSuperExpCutoff.cxx \
                          RadialProfile.cxx \
                          ResponseFunctions.cxx \
                          RoiCuts.cxx \
                          ScData.cxx \
                          ScaleFactor.cxx \
                          SkyDirFunction.cxx \
                          SmoothBrokenPowerLaw.cxx \
                          SmoothDoubleBrokenPowerLaw.cxx \
                          Source.cxx \
                          SourceFactory.cxx \
                          SourceMap.cxx \
                          SourceModel.cxx \
                          SourceModelBuilder.cxx \
                          SpatialMap.cxx \
                          SummedLikelihood.cxx \
                          TiedParameter.cxx \
                          TrapQuad.cxx \
                          WcsMap.cxx \
                          WcsMap2.cxx \
                          WcsMapLibrary.cxx \
                          XmlBuilder.cxx \
                          XmlParser.cxx \
                          dmfit/DMFitFunction.cxx

HD_LIBRARY_SRC_c        = dmfit/dmfit_comm.c \
                          dmfit/dmfit_func.c \
                          dmfit/dmfit_load.c

HD_SUBDIRS              = TsMap \
                          backfile \
                          diffuseResponses \
                          dmfit \
                          expMap \
                          gtaddlivetime \
                          gtexpcube2 \
                          gtfindsrc \
                          gtmodelmap \
                          gtsrcmaps \
                          gtsrcprob \
                          likelihood \
                          makeExposureCube \
                          meanPsf

HD_CXXFLAGS             = -I${HD_SRC}/optimizers \
                          ${HD_STD_CXXFLAGS} \
                          -fpermissive

HD_CLEAN                = ${HD_LIBRARY_SRC_c:.c=.o} \
                          ${HD_LIBRARY_SRC_cxx:.cxx=.o}

HD_CFLAGS               = ${HD_STD_CFLAGS}

HD_INSTALL_LIBRARIES    = ${HD_LIBRARY_ROOT}

include ${HD_STD_MAKEFILE}
