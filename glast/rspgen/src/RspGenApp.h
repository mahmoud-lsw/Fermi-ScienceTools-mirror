/** \file RspGenApp.h
    \brief Declaration for main rspgen application class.
    \author James Peachey, HEASARC
*/
#ifndef rspgen_RspGenApp_h
#define rspgen_RspGenApp_h

#include <string>

#include "evtbin/BinConfig.h"
#include "evtbin/Binner.h"

#include "st_app/StApp.h"

namespace rspgen {

  class IResponse;

  class RspGenApp : public st_app::StApp {
    public:
      RspGenApp();

      virtual ~RspGenApp() throw();

      virtual void run();

      virtual void prompt(st_app::AppParGroup & pars);

      virtual void writeResponse(const st_app::AppParGroup & pars);

      virtual void loadResponses();

    private:
      std::string getDataDir() const;
      evtbin::BinConfig * getConfig(const st_app::AppParGroup & pars);

      evtbin::BinConfig * m_bin_config;
      std::string m_data_dir;
      IResponse * m_response;
  };

}

#endif
