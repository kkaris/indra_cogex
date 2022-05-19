<template>
  <h5 v-if="showTitle">Entity List</h5>
  <div>
    <span
      v-for="([cls, descr], index) in availableClasses"
      class="badge"
      :title="availableClasses.length > 1 ? 'Toggle visibility' : ''"
      @click="toggleHide(cls)"
      :class="`${cls} ${isClsVisible(cls) ? '' : ' opacity-50'}`"
      :key="index"
      >{{ descr }}</span
    >
  </div>
  <hr />
  <div class="text-start">
    <template v-for="(agent, number) in computedList" :key="number">
      <AgentModal :agent-object="agent" />
    </template>
  </div>
</template>

<script>
import AgentModal from "@/components/AgentModal.vue";

export default {
  name: "AgentList.vue",
  components: {
    AgentModal,
  },
  data() {
    return {
      badgeClasses: [
        ["bg-primary", "Gene/Protein"],
        ["bg-secondary", "Small Molecule"],
        ["bg-success", "Biological Process, Disease"],
        ["bg-info text-dark", "Phenotypic Abnormality"],
        ["bg-light text-dark", "Experimental Factor"],
        ["bg-warning text-dark", "ungrounded"],
      ],
      primaryVis: true,
      secondaryVis: true,
      successVis: true,
      infoVis: true,
      lightVis: true,
      warningVis: true,
    };
  },
  props: {
    entities: {
      type: Array,
      required: true,
    },
    showTitle: {
      // Show title
      type: Boolean,
      default: true,
    },
  },
  computed: {
    computedList() {
      let list = [];
      for (let agent of this.entities) {
        if (this.isClsVisible(this.dbRefsToCls(agent.db_refs))) {
          list.push(agent);
        }
      }
      return list;
    },
    availableClasses() {
      let classes = [];
      for (const agent of this.entities) {
        const cls = this.dbRefsToCls(agent.db_refs);
        // Add class if it is in badgeClasses and is not already in classes
        if (
          this.badgeClasses.some((el) => el[0] === cls) &&
          !classes.some((el) => el[0] === cls)
        ) {
          const clsDesc = this.badgeClasses.find((el) => el[0] === cls);
          classes.push(clsDesc);
        }
      }
      return classes;
    },
  },
  methods: {
    toggleHide(cls) {
      // If there is only one class, don't toggle
      if (this.availableClasses.length <= 1) {
        return;
      }
      switch (cls) {
        case "bg-primary":
          this.primaryVis = !this.primaryVis;
          break;
        case "bg-secondary":
          this.secondaryVis = !this.secondaryVis;
          break;
        case "bg-success":
          this.successVis = !this.successVis;
          break;
        case "bg-info text-dark":
          this.infoVis = !this.infoVis;
          break;
        case "bg-light text-dark":
          this.lightVis = !this.lightVis;
          break;
        case "bg-warning text-dark":
          this.warningVis = !this.warningVis;
          break;
      }
    },
    isClsVisible(cls) {
      switch (cls) {
        case "bg-primary":
          return this.primaryVis;
        case "bg-secondary":
          return this.secondaryVis;
        case "bg-success":
          return this.successVis;
        case "bg-info text-dark":
          return this.infoVis;
        case "bg-light text-dark":
          return this.lightVis;
        case "bg-warning text-dark":
          return this.warningVis;
      }
    },
    dbRefsToCls(db_refs) {
      // db_refs is an object with key-value pairs of the form
      // { "namespace": "namespace_id" }
      let cls = "";
      const dbRefKeys = Object.keys(db_refs);
      dbRefKeys.forEach((ns) => {
        // Continue as long as the class is not set or is still default
        if (cls === "warning text-dark" || cls === "") {
          const nsLower = ns.toLowerCase();
          switch (nsLower) {
            case "fplx":
            case "hgnc":
            case "up":
            case "uppro":
            case "mirbase":
              cls = "bg-primary";
              return;
            case "chebi":
              cls = "bg-secondary";
              return;
            case "go":
            case "mesh":
            case "doid":
              cls = "bg-success";
              return;
            case "hp":
              cls = "bg-info text-dark";
              return;
            case "efo":
              cls = "bg-light text-dark";
              return;
            default:
              cls = "warning text-dark";
          }
        }
      });
      return cls;
    },
  },
};
</script>

<style scoped>
.badge {
  cursor: pointer;
}
</style>