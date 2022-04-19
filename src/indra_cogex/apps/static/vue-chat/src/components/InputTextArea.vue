<template>
  <p>Sample textarea</p>
  <textarea :title="title"/>
</template>
<script>
export default {
  name: 'InputTextArea',
  data () {
    return {
      title: 'A title from vue'
    }
  },
}
</script>