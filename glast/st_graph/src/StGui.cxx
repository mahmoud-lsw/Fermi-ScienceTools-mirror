/** \file StAppGui.cxx
    \brief Implementation of StAppGui class.
    \author James Peachey, HEASARC
*/
#include <algorithm>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <utility>

#include "hoops/hoops.h"
//#include "st_app/StApp.h"
//#include "st_app/StAppFactory.h"
//#include "st_app/StEventReceiver.h"

#include "st_graph/Engine.h"
#include "st_graph/IEventReceiver.h"
#include "st_graph/IFrame.h"
#include "st_graph/ITabFolder.h"
#include "st_graph/Placer.h"
#include "st_graph/StGui.h"

#include "st_stream/StreamFormatter.h"

using namespace st_graph;

namespace st_graph {

  ParWidget::ParWidget(Engine & engine, IFrame * parent, hoops::IPar * par, StGui * gui):
    m_engine(engine), m_value_string(), m_gui(gui), m_frame(0), m_label(0), m_value(0), m_open(0),
    m_par(par), m_bool(false), m_stretch(false), m_display(false) {
    if (0 == m_par) throw std::logic_error("ParWidget constructor was passed a null parameter pointer");

    m_frame = m_engine.createComposite(parent, this);

    std::string label = m_par->Name();

    // Get min/max values.
    const std::string & min(m_par->Min());
    const std::string & max(m_par->Max());

    // Handle enumerated list of allowed values.
    if (max.empty() && std::string::npos != min.find("|")) {
      label += " [" + min + "]";
    } else if (!min.empty()) {
      label += " [" + min + ", " + max + "]";
    }

    m_label = m_engine.createLabel(m_frame, this, label);

    // If prompt was supplied, use it to create tool tip.
    const std::string & prompt(m_par->Prompt());
    if (!prompt.empty()) m_label->setToolTipText(prompt);

    // Build width of whole widget from constituent widths.
    // Width = width of label + ...
    long width = m_label->getWidth();

    if (std::string::npos != m_par->Type().find("b")) {
      m_bool = true;
      // If boolean parameter, use a checkbox.
      m_value = m_engine.createButton(m_frame, this, "check", "");

      // Set button in state consistent with parameter value.
      bool state = *m_par;
      if (state) {
        m_value_string = "down";
      } else {
        m_value_string = "up";
      }
      m_value->setState(m_value_string);
    } else {
      // Get value from parameter.
      m_value_string = m_par->Value();

      // For all other non-boolean parameters, use a text edit.
      m_value = m_engine.createTextEntry(m_frame, this, m_value_string);

      // Set width to the standard width for this parameter type.
      m_value->setWidth(entryWidth(m_par));

      // Special cases.
      if (std::string::npos != m_par->Type().find("f")) {
        // File types have additional option of a file dialog box, activated by an "Open" button.
        m_open = m_engine.createButton(m_frame, this, "text", "...");
        // Adjust width to include button.
        width += m_open->getWidth() + 3;
        m_stretch = true;
      } else if (std::string::npos != m_par->Type().find("s")) {
        m_stretch = true;
      }
    }
    // Adjust width to include m_value widget.
    width += m_value->getWidth() + 3;

    // Make certain widget will not shrink smaller than the label size.
    m_frame->setMinimumWidth(m_label->getWidth());

    m_frame->setWidth(width);
    m_frame->setHeight(std::max(m_label->getHeight(), m_value->getHeight()));

    // By default, "hidden" parameters are not displayed to start with.
    m_display = (std::string::npos == m_par->Mode().find("h"));
  }

  ParWidget::~ParWidget() { delete m_frame; }

  void ParWidget::layout(IFrame * f) {
    if (f != m_frame) return;

    Center(m_value).below(Center(m_frame));
    Center(m_label).below(Center(m_frame));
    if (0 != m_open) Center(m_open).below(Center(m_frame));

    LeftEdge(m_label).rightOf(LeftEdge(m_frame));
    LeftEdge(m_value).rightOf(RightEdge(m_label), 3);

    if (0 == m_open) {
      if (m_stretch) RightEdge(m_value).stretchTo(RightEdge(m_frame));
    } else {
      RightEdge(m_open).leftOf(RightEdge(m_frame));
      if (m_stretch) RightEdge(m_value).stretchTo(RightEdge(m_open), -3);
    }

    if (!m_display) {
      m_frame->unDisplay();
      m_label->unDisplay();
      m_value->unDisplay();
      if (0 != m_open) m_open->unDisplay();
    }
  }

  void ParWidget::clicked(IFrame * f) {
    std::string state;
    if (m_open == f && m_open != 0) {
      if (std::string::npos != m_par->Type().find("w"))
        state = m_engine.fileDialog(m_frame, m_value_string, "save");
      else
        state = m_engine.fileDialog(m_frame, m_value_string, "open");
    } else {
      state = m_value->getState();
    }
    m_gui->synchronizeWidgets(getName(), state);
  }

  void ParWidget::modified(IFrame *, const std::string & text) {
    m_gui->synchronizeWidgets(getName(), text);
  }

  ParWidget::operator IFrame * () { return getFrame(); }

  IFrame * ParWidget::getFrame () { return m_frame; }

  IFrame * ParWidget::getLabel() { return m_label; }

  const std::string & ParWidget::getName() const { return m_par->Name(); }

  std::string ParWidget::getValue() const {
    if (m_bool) {
      if (m_value_string == "down") return "true";
      return "false";
    }
    return m_value_string;
  }

  void ParWidget::setValue(const std::string & value_string) {
    m_value_string = value_string;
    m_value->setState(value_string);
  }

  void ParWidget::display(bool disp_flag) {
    m_display = disp_flag;
    if (m_display) {
      if (0 != m_open) m_open->display();
      m_value->display();
      m_label->display();
      m_frame->display();
    } else {
      layout(m_frame);
    }
  }

  long ParWidget::entryWidth(hoops::IPar * par) const {
    // The first time this is called, create a temporary gui with text entries corresponding to the sizes of parameters.
    static std::map<std::string, long> s_width;
    if (s_width.empty()) {
      std::auto_ptr<IFrame> mf(m_engine.createMainFrame(0, 100, 100, "sizer"));
      std::auto_ptr<IFrame> bool_pw(m_engine.createTextEntry(mf.get(), 0, "false"));
      std::auto_ptr<IFrame> int_pw(m_engine.createTextEntry(mf.get(), 0, "+1234567890"));
      std::auto_ptr<IFrame> float_pw(m_engine.createTextEntry(mf.get(), 0, "1.2345678901234E+123"));
      std::auto_ptr<IFrame> string_pw(m_engine.createTextEntry(mf.get(), 0, "1234567890123456789012345678901234567890"));
      //std::auto_ptr<IFrame> string_pw(m_engine.createTextEntry(mf.get(), 0, "123456789012345678901234"));
      // Store sizes of text entry boxes for each parameter type.
      s_width.insert(std::make_pair(std::string("b"), bool_pw->getWidth()));
      s_width.insert(std::make_pair(std::string("i"), int_pw->getWidth()));
      s_width.insert(std::make_pair(std::string("r"), float_pw->getWidth()));
      s_width.insert(std::make_pair(std::string("f"), string_pw->getWidth()));
      s_width.insert(std::make_pair(std::string("s"), string_pw->getWidth()));
    }

    long width = 10;
    for (std::map<std::string, long>::iterator itor = s_width.begin(); itor != s_width.end(); ++itor) {
      if (std::string::npos != par->Type().find(itor->first)) {
        width = itor->second;
        break;
      }
    }
    return width;
  }

  hoops::IParGroup & operator <<(hoops::IParGroup & group, const ParWidget & par) {
    group[par.getName()] = par.getValue();
    return group;
  }

  //StEventReceiver::StEventReceiver(st_graph::Engine & engine, hoops::IParGroup & par_group, StEventReceiver * app):
  StGui::StGui(Engine & engine, const hoops::IParGroup & par_group):
    m_os("StGui", "StGui", 2), m_engine(engine), m_par_widget(), m_tab_folder(), m_parent(), m_plot_title(),
    m_par_group(par_group.Clone()), m_main(0), m_group_frame(0), m_run(0), m_cancel(0), m_show_advanced(0),
    m_plot_frame(0), m_widest(0), m_tab_height(0), m_plot_enabled(false) {
    try {
      m_plot_enabled = (*m_par_group)["plot"];
    } catch (const std::exception &) {
      // Ignore this exception.
    }
    try {
      m_plot_title = (*m_par_group)["title"].Value();
    } catch (const std::exception &) {
      // Ignore this exception.
    }

    // Prevent most windows from exiting this application.
    engine.setDefaultExitOnClose(false);
  }

  StGui::~StGui() {
    for (ParWidgetCont::reverse_iterator itor = m_par_widget.rbegin(); itor != m_par_widget.rend(); ++itor)
      delete *itor->second;
    delete m_main;
    delete m_par_group;
  }

  void StGui::clicked(IFrame * f) {
    hoops::IParGroup & pars(*m_par_group);
    if (f == m_run) {

      try {
        // Get parameter values which are associated with the state of a tab folder.
        for (TabFolderCont::iterator itor = m_tab_folder.begin(); itor != m_tab_folder.end(); ++itor) {
//          pars[itor->first] = itor->second->getSelected();
        }

        // Get parameter values which are associated with parameter widgets.
        for (ParWidgetCont::iterator itor = m_par_widget.begin(); itor != m_par_widget.end(); ++itor) {
          pars << *itor->second;
        }

//        pars.Save();
      } catch (const std::exception & x) {
        m_os.err() << "Problem with parameter: " << x.what() << std::endl;
        return;
      }

      runApp();
#if 0
      try {
        int chatter = pars["chatter"];
        IStAppFactory::instance().setMaximumChatter(chatter);
      } catch (const std::exception &) {
        // Ignore
      }

      try {
        bool debug = pars["debug"];
        IStAppFactory::instance().setDebugMode(debug);
      } catch (const std::exception &) {
        // Ignore
      }

      try {
        m_app->run();
      } catch (const std::exception & x) {
        m_os.err() << "Running the application failed: " << std::endl << x.what() << std::endl;
      }
#endif
    } else if (f == m_cancel) {
      m_engine.stop();
    } else if (f == m_show_advanced) {
      // Show/hide advanced parameters.
      for (ParWidgetCont::iterator itor = m_par_widget.begin(); itor != m_par_widget.end(); ++itor) {
        // If parameter is not "hidden" do not tamper with its visibility.
        if (std::string::npos == pars[itor->first].Mode().find("h")) continue;
        const std::string & state(m_show_advanced->getState());
        if (state == "up") itor->second->display(false);
        else if (state == "down") itor->second->display(true);
      }
      layout(m_group_frame);
    }
  }

  void StGui::closeWindow(IFrame * f) {
    if (f == m_main) m_engine.stop();
  }

  void StGui::layout(IFrame * f) {
    if (f == m_main) {
      // Stack buttons horizontally at the top of the frame.
      LeftEdge(m_run).rightOf(LeftEdge(m_main), 6);
      LeftEdge(m_cancel).rightOf(RightEdge(m_run));
      LeftEdge(m_show_advanced).rightOf(LeftEdge(m_run));

      TopEdge(m_run).below(TopEdge(m_main), 6);
      TopEdge(m_cancel).below(TopEdge(m_main), 6);
      TopEdge(m_show_advanced).below(BottomEdge(m_cancel), 6);

      // Size the group frame so that it sits nicely below the buttons.
      TopEdge(m_group_frame).below(BottomEdge(m_show_advanced), 6);
//      BottomEdge(m_group_frame).stretchTo(BottomEdge(m_main), -6);
      LeftEdge(m_group_frame).rightOf(LeftEdge(m_main), 6);
      RightEdge(m_group_frame).stretchTo(RightEdge(m_main), -6);
//      RightEdge(m_group_frame).stretchTo(RightEdge(m_widest->getFrame()));

      // Size the plot so it sits nicely to the right of the group frame, and maintains constant aspect ratio.
      if (0 != m_plot_frame) {
        TopEdge(m_plot_frame).below(BottomEdge(m_group_frame), 6);
//        TopEdge(m_plot_frame).below(BottomEdge(m_show_advanced), 6);
        LeftEdge(m_plot_frame).rightOf(LeftEdge(m_main), 6);
//        LeftEdge(m_plot_frame).rightOf(RightEdge(m_group_frame), 6);
        RightEdge(m_plot_frame).stretchTo(RightEdge(m_main), -6);
        BottomEdge(m_plot_frame).stretchTo(BottomEdge(m_main), -6);
//        m_plot_frame->setHeight(m_plot_frame->getWidth() / 2);
      }
    
      // Layout tab folders.
      for (TabFolderCont::iterator tab_itor = m_tab_folder.begin(); tab_itor != m_tab_folder.end(); ++tab_itor) {
        // Starting height of tab folders is 0.
        long height = 0;
        // Fill a container with this folder's tabs.
        std::map<std::string, IFrame *> tab_cont;
        tab_itor->second->getTabCont(tab_cont);

        // Layout widgets on each tab.
        for (std::map<std::string, IFrame *>::iterator itor = tab_cont.begin(); itor != tab_cont.end(); ++itor) {
          layout(itor->second);
          height = height > itor->second->getHeight() ? height : itor->second->getHeight();
        }
        // Set overall height of the folder to the computed height + the height of the tabs.
        tab_itor->second->getFrame()->setHeight(height + m_tab_height + 10);
      }
    } else {
      std::list<IFrame *> subframes;
      f->getSubframes(subframes);

      std::list<IFrame *>::iterator itor = subframes.begin();
      if (itor != subframes.end()) {
        IFrame * previous = *itor;
        TopEdge(*itor).below(TopEdge(f), 22);
        LeftEdge(*itor).rightOf(LeftEdge(f), 10);
        RightEdge(*itor).stretchTo(RightEdge(f), -10);

        for (++itor; itor != subframes.end(); ++itor) {
          TopEdge(*itor).below(BottomEdge(previous), 6);
          LeftEdge(*itor).rightOf(LeftEdge(f), 10);
          RightEdge(*itor).stretchTo(RightEdge(f), -10);
          previous = *itor;
        }
        BottomEdge(f).stretchTo(BottomEdge(previous), 10);
        BottomEdge(m_group_frame).stretchTo(BottomEdge(f));
      }
    }
  }

  void StGui::run() {
    // Set up standard Gui main window.
    createMainFrame();

    hoops::IParGroup & pars(*m_par_group);
    for (hoops::GenParItor itor = pars.begin(); itor != pars.end(); ++itor) {
      const std::string & par_name((*itor)->Name());

      // Changing from GUI to command line mode is not permitted. Also, mode is irrelevant.
      // Skip blank lines as well.
      if (par_name == "gui" || par_name == "mode") continue;
      else if (0 == par_name.size()) continue;

      // Get range of parameter.
      std::list<std::string> par_range;
//      bool enumerated_range = parseRange(*itor, par_range);

      // Get the appropriate parent frames for this parameter.
      std::list<IFrame *> parent;
      getParent(*itor, parent);

      // Loop over parents.
      for (std::list<IFrame *>::iterator parent_itor = parent.begin(); parent_itor != parent.end(); ++parent_itor) {
        // If parameter is a switch with an enumerated range, make it a tab-folder.
        //if (m_par_group.isSwitch(par_name) && enumerated_range) {
        if (false) {
          ITabFolder * tf = m_engine.createTabFolder(*parent_itor, this);
          for (std::list<std::string>::iterator enum_itor = par_range.begin(); enum_itor != par_range.end(); ++enum_itor) {
            IFrame * frame = tf->addTab(*enum_itor);
            // Use the current parameter value to select the correct tab.
            if ((*itor)->Value() == *enum_itor) tf->select(frame);
          }
          tf->getFrame()->setNaturalSize();
          m_tab_height = tf->getFrame()->getHeight();
          m_tab_folder.insert(std::make_pair(par_name, tf));

          // Record parent of this widget.
          m_parent.insert(std::make_pair(tf->getFrame(), *parent_itor));
        } else {
          // Create standard widget representing each parameter.
          ParWidget * widget = createParWidget(*itor, *parent_itor);

          // Store widget in container.
          m_par_widget.insert(std::make_pair(par_name, widget));

          // Keep track of the widget with the widest label.
          if (0 == m_widest || widget->getLabel()->getWidth() > m_widest->getLabel()->getWidth()) m_widest = widget;

          // Record parent of this widget.
          m_parent.insert(std::make_pair(widget->getFrame(), *parent_itor));
        }
      }
    }

    if (0 != m_widest) {
      m_group_frame->setMinimumWidth(m_widest->getFrame()->getWidth() + 12);
      for (ParWidgetCont::iterator itor = m_par_widget.begin(); itor != m_par_widget.end(); ++itor) {
        itor->second->getLabel()->setWidth(m_widest->getLabel()->getWidth());
      }
    }
    if (0 != m_plot_frame) m_plot_frame->setMinimumWidth(100);

    m_engine.run();
  }

  void StGui::createMainFrame() {
    // Use the name and version of the tool as a label for the GUI window.
#if 0
    std::string label(m_app->getName());
    if (!label.empty()) label += " ";
    const std::string & version(m_app->getVersion());
    if (!version.empty()) label += "version " + version;
#endif
    std::string label("StGui");

    m_main = m_engine.createMainFrame(this, 650, 600, label);
    m_group_frame = m_engine.createGroupFrame(m_main, this, "Parameters");
    m_run = m_engine.createButton(m_main, this, "text", "Run");
    m_cancel = m_engine.createButton(m_main, this, "text", "Cancel");
    m_show_advanced = m_engine.createButton(m_main, this, "check", "Show Advanced Parameters");
    if (m_plot_enabled) m_plot_frame = m_engine.createPlotFrame(m_main, m_plot_title, 638, 319);

    // Set up some tool tips.
    m_run->setToolTipText("Run the application from inside the GUI");
    m_cancel->setToolTipText("Exit the application and GUI");
    m_show_advanced->setToolTipText("Display advanced (\"hidden\") parameters");

    // Disable prompting.
    //hoops::IParGroup & pars(m_par_group);
    //pars.setPromptMode(false);
  }

  ParWidget * StGui::createParWidget(hoops::IPar * par, IFrame * parent) {
    return new ParWidget(m_engine, parent, par, this);
  }

  void StGui::synchronizeWidgets(const std::string & par_name, const std::string & value) {
    std::pair<ParWidgetCont::iterator, ParWidgetCont::iterator> range = m_par_widget.equal_range(par_name);
    for (ParWidgetCont::iterator itor = range.first; itor != range.second; ++itor) {
      itor->second->setValue(value);
    }
  }

  void StGui::enablePlotFrame(const std::string & title) {
    m_plot_enabled = true;
    m_plot_title = title;
  }

  IFrame * StGui::getPlotFrame() { return m_plot_frame; }

  const IFrame * StGui::getPlotFrame() const { return m_plot_frame; }

  bool StGui::parseRange(const hoops::IPar * par, std::list<std::string> & range) {
    range.clear();

    const std::string & par_min(par->Min());
    const std::string & par_max(par->Max());

    // Check whether min/max is really min/max or defines a set of enumerated possible values.
    bool enumerated_range = (par_max.empty() && std::string::npos != par_min.find("|"));

    if (!enumerated_range) {
      // min/max are simply min/max.
      range.push_back(par_min);
      range.push_back(par_max);
    } else {
      // Parse enumerated range.
      std::string::const_iterator begin = par_min.begin();
      std::string::const_iterator end = par_min.end();
      while (begin != end) {
        // Skip leading whitespace.
        while (begin != end && isspace(*begin)) ++begin;

        // Move to end of token (end of string or | or whitespace)
        std::string::const_iterator itor = begin;
        for (; itor != end && '|' != *itor && !isspace(*itor); ++itor) {}

        // Save this token in output enumerated range container.
        if (begin != itor) {
          range.push_back(std::string(begin, itor));
        }

        // Skip trailing whitespace.
        for (begin = itor; begin != end && isspace(*begin); ++begin) {}

        // Skip |s.
        while (begin != end && '|' == *begin) ++begin;
      }
    }

    return enumerated_range;
  }

  //void StEventReceiver::getParent(const hoops::IPar * par, std::list<st_graph::IFrame *> & parent) {
  void StGui::getParent(const hoops::IPar *, std::list<IFrame *> & parent) {
//    const std::string & name(par->Name());
    // Clear out previous container of frames.
    parent.clear();

    // Get cases on which this parameter depends.
#if 0
    AppParGroup::CaseList case_cont;
    m_par_group.getCase(name, case_cont);

    // Loop over all cases.
    for (AppParGroup::CaseList::iterator itor = case_cont.begin(); itor != case_cont.end(); ++itor) {
      // itor->first == name of switch.
      // itor->second == value of switch.
      // See if switch is displayed in one or more tab folders.
      std::pair<TabFolderCont::iterator, TabFolderCont::iterator> range = m_tab_folder.equal_range(itor->first);
      for (TabFolderCont::iterator tab_itor = range.first; tab_itor != range.second; ++tab_itor) {
        // Find the tab corresponding to the parameter value given by the second part of the case.
        IFrame * frame = tab_itor->second->getTab(itor->second);
        if (0 != frame) parent.push_back(frame);
      }
    }
#endif

    if (parent.empty()) parent.push_back(m_group_frame);
  }

}
