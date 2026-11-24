/* DMTCP hooks: copy/overwrite only regular, non-executable study files.
 * Keep global --checkpoint-open-files and --allow-file-overwrite DISABLED.
 * Executables, input data and dependencies are reused in place, never replaced.
 */
#define _GNU_SOURCE
#include <dmtcp.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static int owned_file(const char *path)
{
    const char *root = getenv("NEPAL_DMTCP_FILE_ROOT");
    char resolved[PATH_MAX];
    struct stat info;
    if (!root || root[0] != '/' || !path || !realpath(path, resolved)) return 0;
    size_t n = strlen(root);
    if (n < 2 || strncmp(resolved, root, n) || resolved[n] != '/') return 0;
    if (stat(resolved, &info) || !S_ISREG(info.st_mode)) return 0;
    if (info.st_mode & (S_IXUSR | S_IXGRP | S_IXOTH)) return 0;
    return 1;
}

int dmtcp_must_ckpt_file(const char *path) { return owned_file(path); }
int dmtcp_must_overwrite_file(const char *path) { return owned_file(path); }

static void event_hook(DmtcpEvent_t event, DmtcpEventData_t *data)
{
    (void)event;
    (void)data;
}

static DmtcpPluginDescriptor_t policy = {
    DMTCP_PLUGIN_API_VERSION, DMTCP_PACKAGE_VERSION,
    "nepal-study-files", "Nepal simulation project", "",
    "Checkpoint only study-owned data; preserve installed executables", event_hook
};
DMTCP_DECL_PLUGIN(policy);
