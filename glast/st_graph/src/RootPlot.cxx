/** \file RootPlot.cxx
    \brief Implementation of RootPlot class.
    \author James Peachey, HEASARC/GSSC
*/
#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "RootPlot.h"
#include "RootPlotFrame.h"

#include "st_graph/IFrame.h"
#include "st_graph/Sequence.h"

namespace st_graph {

  RootPlot::RootPlot(IFrame * parent, const std::string & style, const ISequence & x, const ISequence & y, bool delete_parent):
    m_seq_cont(0), m_label(), m_style(), m_line_style("solid"), m_curve_type("line"), m_line_color(Color::eBlack), m_dimensionality(2),
    m_parent(0), m_z_data(0), m_delete_parent(delete_parent) {
    // Get the parent multi frame so that the plot can be added with desired style.
    m_parent = dynamic_cast<RootPlotFrame *>(parent);
    if (0 == m_parent) throw std::logic_error("RootPlot constructor: parent must be a valid RootPlotFrame");

    // Sanity check.
    if (x.size() != y.size()) throw std::logic_error("RootPlot constructor: x and y sequence do not have same size");

    // Add this plot to parent's container of plots, allowing for auto-delete.
    m_parent->addPlot(this);

    setStyle(style);

    m_seq_cont.push_back(x.clone());
    m_seq_cont.push_back(y.clone());
  }

  RootPlot::RootPlot(IFrame * parent, const std::string & style, const ISequence & x, const ISequence & y,
    const std::vector<std::vector<double> > & z, bool delete_parent): m_seq_cont(0), m_label(), m_style(),
    m_line_style("solid"), m_curve_type("line"), m_line_color(Color::eBlack), m_dimensionality(3), m_parent(0),
    m_z_data(&z), m_delete_parent(delete_parent) {
    // Get the parent multi frame so that the plot can be added with desired style.
    m_parent = dynamic_cast<RootPlotFrame *>(parent);
    if (0 == m_parent) throw std::logic_error("RootPlot constructor: parent must be a valid RootPlotFrame");

    // Sanity check.
    if (x.size() != z.size())
      throw std::logic_error("RootPlot constructor: x sequence and first data dimension do not have same size");
    if (y.size() != z.begin()->size())
      throw std::logic_error("RootPlot constructor: y sequence and second data dimension do not have same size");

    // Add this plot to parent's container of plots, allowing for auto-delete.
    m_parent->addPlot(this);

    setStyle(style);

    m_seq_cont.push_back(x.clone());
    m_seq_cont.push_back(y.clone());
  }

  RootPlot::~RootPlot() {
    // Note: This appears more complicated than necessary, but be careful changing it. Under some circumstances,
    // a RootPlot needs to delete its parent, but the parent will always attempt to delete the RootPlot in the
    // process. Thus it is important to ensure the child is detached at the right times to prevent deleting
    // the parent and/or the child twice.

    // Save pointer to parent in case delete parent was selected.
    IFrame * parent = m_parent;

    // Break all links between this and its parent. After this, m_parent == 0.
    if (0 != m_parent) m_parent->removePlot(this);

    // Delete sequences.
    for (std::vector<const ISequence *>::reverse_iterator itor = m_seq_cont.rbegin(); itor != m_seq_cont.rend(); ++itor)
      delete (*itor);

    // In special case where the plot owns its parent frame, delete that as well.
    if (m_delete_parent) delete parent;
  }

  const std::vector<const ISequence *> RootPlot::getSequences() const { return m_seq_cont; }

  const std::vector<std::vector<double> > & RootPlot::getZData() const {
    if (0 == m_z_data) throw std::logic_error("RootPlot::getZData() called for a plot which has null Z data");
    return *m_z_data;
  }

  std::vector<Axis> & RootPlot::getAxes() { return m_parent->getAxes(); }

  const std::vector<Axis> & RootPlot::getAxes() const { return m_parent->getAxes(); }

  void RootPlot::addMarker(double x, double y, const std::string & text, int color) {
    m_label.push_back(Marker(x, y, text, color));
    m_parent->addMarker(m_label.back());
  }

  void RootPlot::getMarkers(std::vector<Marker> & labels) const {
    labels.assign(m_label.begin(), m_label.end());
  }

  std::vector<Marker> & RootPlot::getMarkers() { return m_label; }

  int RootPlot::getLineColor() const { return m_line_color; }

  void RootPlot::setLineColor(int color) { m_line_color = color; }

  std::string RootPlot::getLineStyle() const { return m_line_style; }

  void RootPlot::setLineStyle(const std::string & style) { m_line_style = style; }

  std::string RootPlot::getCurveType() const { return m_curve_type; }

  void RootPlot::setCurveType(const std::string & type) { m_curve_type = type; }

  const std::string & RootPlot::getStyle() const { return m_style; }

  void RootPlot::setStyle(const std::string & style) {
    m_style = style;

    // Convert style string to lower case.
    for (std::string::iterator itor = m_style.begin(); itor != m_style.end(); ++itor) *itor = tolower(*itor);

    if (std::string::npos != m_style.find("h")) {
      m_style = "hist";
    } else if (std::string::npos != m_style.find("l")) {
      m_style = "lego";
    } else if (std::string::npos != m_style.find("sc")) {
      m_style = "scat";
    } else if (std::string::npos != m_style.find("su")) {
      m_style = "surf";
    } else {
      throw std::logic_error("RootPlot::setStyle: unknown plot style \"" + style + "\"");
    }

  }

  void RootPlot::setParent(RootPlotFrame * parent) { m_parent = parent; }

}
