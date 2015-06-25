// $Id: ScheduledEvent.cxx,v 1.1.1.1 1999/12/20 22:27:42 burnett Exp $


#include "facilities/ScheduledEvent.h"
#include "facilities/Scheduler.h"

#include <typeinfo>
void ScheduledEvent::schedule(double t, ScheduledEvent* next)
{
    Scheduler::instance()->schedule(t, next);
}

std::string ScheduledEvent::name()const
{
    const std::type_info& t = typeid(*this);
    return t.name();
}
