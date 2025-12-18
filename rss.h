//
// Created by Max Shikunov on 18/12/2025.
//

#ifndef UNTITLED3_RSS_H
#define UNTITLED3_RSS_H
#include <mach/mach.h>
#include <cstdint>


uint64_t current_rss_bytes() {
    mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;

    if (task_info(mach_task_self(),
                  MACH_TASK_BASIC_INFO,
                  (task_info_t)&info,
                  &count) != KERN_SUCCESS)
        return 0;

    return info.resident_size;
}

#endif //UNTITLED3_RSS_H
